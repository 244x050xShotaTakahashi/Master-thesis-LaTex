import argparse
import csv
import shutil
from pathlib import Path
from typing import List, Optional, Union

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Circle

# 入力ファイルのデフォルトを変更 (CSV)
DATA_FILE_DEFAULT = Path(__file__).resolve().parent.parent  / "cutoff_0.01_bimodal" / "job_6137618_task_7" / "particles.csv"
WALLS_FILE_DEFAULT = Path(__file__).resolve().parent.parent / "DEM-valid" / "inputs" / "walls.dat"

def read_walls_data(filename: Union[str, Path] = WALLS_FILE_DEFAULT):
    """walls.dat を読み込み、斜面壁のリストを返す。
    
    各行は x_start z_start x_end z_end の形式。
    コメント行（#や!で始まる）や空行はスキップされる。
    """
    walls = []
    try:
        with open(filename, "r") as fh:
            for line in fh:
                line = line.strip()
                # コメント行や空行をスキップ
                if not line or line.startswith('#') or line.startswith('!'):
                    continue
                
                # コメント部分を削除
                if '#' in line:
                    line = line.split('#')[0].strip()
                if '!' in line:
                    line = line.split('!')[0].strip()
                
                try:
                    values = [float(x) for x in line.split()]
                    if len(values) >= 4:
                        walls.append({
                            'x_start': values[0],
                            'z_start': values[1],
                            'x_end': values[2],
                            'z_end': values[3]
                        })
                except ValueError:
                    continue  # 解析できない行はスキップ
    except FileNotFoundError:
        pass
    
    return walls

def read_simulation_data(
    filename: Union[str, Path] = DATA_FILE_DEFAULT,
    frame_step: int = 1,
    max_frames: Optional[int] = None,
):
    """particles.csv を読み込み、各タイムステップのデータを辞書のリストとして返す。
    
    CSVヘッダー: step,time,id,x,z,radius,vx,vz,omega,angle,mass,charge
    """

    if frame_step < 1:
        raise ValueError("frame_step must be >= 1")
    if max_frames is not None and max_frames < 1:
        raise ValueError("max_frames must be >= 1 when provided")

    frames_data = []
    print("データファイルを読み込み中...")
    
    # データをステップごとにグループ化するための一時辞書
    # key: step, value: {time, particles: [...]}
    temp_frames = {}
    
    try:
        with open(filename, "r", newline='') as fh:
            reader = csv.DictReader(fh)
            
            # ヘッダーのチェック (最低限必要なカラムがあるか)
            required_columns = {'step', 'time', 'x', 'z', 'radius'}
            if not required_columns.issubset(reader.fieldnames):
                print(f"エラー: CSVファイルに必要な列が含まれていません。必要: {required_columns}")
                return None

            for row in reader:
                try:
                    step = int(row['step'])
                    
                    # ステップ間引き (読み込み時に判定してメモリ節約)
                    # 注意: CSVは行ごとなので、新しいステップが現れたタイミングで判定する等の工夫が必要だが、
                    # ここでは単純に全データを舐めて、後で間引くか、step番号で判定する。
                    # データ量が膨大な場合、step % frame_step != 0 なら continue するのが効率的。
                    
                    # CSVのstepは1から始まると仮定。frame_stepごとのデータを取得。
                    # step自体が飛び飛びの場合もあるので、出現順序（frame_count）での間引きは難しい。
                    # ここでは step 番号に基づいてフィルタリングする。
                    if (step - 1) % frame_step != 0:
                        continue
                        
                    time_val = float(row['time'])
                    
                    if step not in temp_frames:
                        if max_frames is not None and len(temp_frames) >= max_frames:
                            # 最大フレーム数に達したら、新しいステップの追加はやめる
                            # ただし、現在のステップが既存なら追加し続ける必要がある
                            # ここでは「読み込み済みステップ数」で判定
                            # (CSVがstep順に並んでいると仮定)
                            break
                        
                        temp_frames[step] = {
                            'time': time_val,
                            'particles': [],
                            # width/height はCSVに含まれていないため、後で計算または固定値を使用
                            'container_width': 0.0, 
                            'container_height': 0.0 
                        }
                    
                    # 粒子データの追加
                    p_data = {
                        "x": float(row['x']),
                        "z": float(row['z']),
                        "r": float(row['radius']),
                        "rotation_angle": float(row['angle']) if 'angle' in row else 0.0,
                    }
                    temp_frames[step]['particles'].append(p_data)

                except ValueError as e:
                    continue # 数値変換エラー等はスキップ

    except FileNotFoundError:
        print(f"エラー: ファイル '{filename}' が見つかりません。")
        return None
    except Exception as e:
        print(f"ファイル読み込み中にエラーが発生しました: {e}")
        return None
    
    # 辞書からリストへ変換（ステップ順にソート）
    sorted_steps = sorted(temp_frames.keys())
    
    for step in sorted_steps:
        frame = temp_frames[step]
        particles = frame['particles']
        num_particles = len(particles)
        
        # コンテナサイズの推定（CSVに含まれていないため、粒子の最大位置から推定）
        # 必要に応じてコマンドライン引数で渡すように変更することも可能
        max_x = max((p['x'] + p['r'] for p in particles), default=1.0)
        max_z = max((p['z'] + p['r'] for p in particles), default=1.0)
        
        # 余裕を持たせる
        container_width = max_x * 1.1
        container_height = max_z * 1.1
        
        frames_data.append({
            "time": frame['time'],
            "num_particles": num_particles,
            "container_width": container_width,   # 推定値
            "container_height": container_height, # 推定値
            "particles": particles,
        })

    print(f"\r  {len(frames_data)} フレームの読み込み完了!        ")
    return frames_data

def animate(frames_data, output_filename="pem_animation.mp4", walls_data=None, fps: int = 10):
    if not frames_data:
        print("アニメーションするデータがありません。")
        return

    # フォント設定を無効化して英語のみ使用
    import matplotlib
    matplotlib.rcParams['font.family'] = 'DejaVu Sans'

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel("X coordinate")
    ax.set_ylabel("Z coordinate")
    fig.suptitle("PEM Validation", fontsize=12)

    if walls_data is None:
        walls_data = []

    # 軸範囲（粒子 + 壁の位置を考慮）
    particle_min_x = None
    particle_max_x = None
    particle_min_z = None
    particle_max_z = None

    for frame_data in frames_data:
        for p in frame_data['particles']:
            px_min = p['x'] - p['r']
            px_max = p['x'] + p['r']
            pz_min = p['z'] - p['r']
            pz_max = p['z'] + p['r']

            particle_min_x = px_min if particle_min_x is None else min(particle_min_x, px_min)
            particle_max_x = px_max if particle_max_x is None else max(particle_max_x, px_max)
            particle_min_z = pz_min if particle_min_z is None else min(particle_min_z, pz_min)
            particle_max_z = pz_max if particle_max_z is None else max(particle_max_z, pz_max)

    if particle_min_x is None:
        particle_min_x = 0.0
    if particle_max_x is None:
        particle_max_x = 1.0
    if particle_min_z is None:
        particle_min_z = 0.0
    if particle_max_z is None:
        particle_max_z = 1.0

    wall_min_x = None
    wall_max_x = None
    wall_min_z = None
    wall_max_z = None
    for wall in walls_data or []:
        xmin = min(wall['x_start'], wall['x_end'])
        xmax = max(wall['x_start'], wall['x_end'])
        zmin = min(wall['z_start'], wall['z_end'])
        zmax = max(wall['z_start'], wall['z_end'])
        wall_min_x = xmin if wall_min_x is None else min(wall_min_x, xmin)
        wall_max_x = xmax if wall_max_x is None else max(wall_max_x, xmax)
        wall_min_z = zmin if wall_min_z is None else min(wall_min_z, zmin)
        wall_max_z = zmax if wall_max_z is None else max(wall_max_z, zmax)

    min_x = particle_min_x if wall_min_x is None else min(particle_min_x, wall_min_x)
    max_x = particle_max_x if wall_max_x is None else max(particle_max_x, wall_max_x)
    min_z = particle_min_z if wall_min_z is None else min(particle_min_z, wall_min_z)
    max_z = particle_max_z if wall_max_z is None else max(particle_max_z, wall_max_z)

    width = max_x - min_x
    height = max_z - min_z
    if width == 0:
        width = 1.0
    if height == 0:
        height = 1.0

    x_margin = 0.1 * width
    y_margin = 0.1 * height

    x_plot_min = min_x - x_margin
    x_plot_max = max_x + x_margin
    y_plot_min = min_z - y_margin
    y_plot_max = max_z + y_margin

    ax.set_xlim(x_plot_min, x_plot_max)
    ax.set_ylim(y_plot_min, y_plot_max)

    time_text_obj = ax.text(0.05, 0.95, '', transform=ax.transAxes, ha="left", va="top", fontsize=10)

    # 壁（簡易表示: 床のみ、あるいは矩形）
    # CSVには壁情報がないため、床(z=0)と左右の境界線を表示する程度にする
    floor_line = Line2D([-100, 100], [0, 0], color='black', lw=2)
    floor_line.set_animated(True)
    ax.add_line(floor_line)
    
    wall_lines: List[Line2D] = []
    if walls_data:
        for idx, wall in enumerate(walls_data):
            line = Line2D(
                [wall['x_start'], wall['x_end']],
                [wall['z_start'], wall['z_end']],
                color='black',
                lw=3.0,
                alpha=0.9,
                zorder=5,
            )
            ax.add_line(line)
            line.set_animated(True)
            wall_lines.append(line)

    vertical_wall_lines: List[Line2D] = []
    if particle_min_x is not None and particle_max_x is not None:
        for x_pos in (particle_min_x, particle_max_x):
            line = Line2D(
                [x_pos, x_pos],
                [y_plot_min, y_plot_max],
                color='black',
                lw=2.5,
                linestyle='-',
                alpha=0.85,
                zorder=4,
            )
            ax.add_line(line)
            line.set_animated(True)
            vertical_wall_lines.append(line)
            
    max_particles = max((frame['num_particles'] for frame in frames_data), default=0)
    particle_patches: List[Circle] = []
    rotation_lines: List[Line2D] = []
    
    for _ in range(max_particles):
        circle = Circle((0.0, 0.0), 0.0, facecolor='white', edgecolor='black', linewidth=1.5, alpha=0.9)
        circle.set_visible(False)
        ax.add_patch(circle)
        particle_patches.append(circle)

        rotation_line = Line2D([], [], color='black', linewidth=1, alpha=0.8)
        rotation_line.set_visible(False)
        ax.add_line(rotation_line)
        rotation_lines.append(rotation_line)

    artists_to_update: List[object] = [time_text_obj, floor_line]
    artists_to_update.extend(wall_lines)
    artists_to_update.extend(vertical_wall_lines)
    artists_to_update.extend(particle_patches)
    artists_to_update.extend(rotation_lines)

    frame_counter = {'count': 0}
    total_frames = len(frames_data)

    def update_frame(frame_idx):
        frame_counter['count'] += 1
        if frame_counter['count'] % 10 == 0 or frame_counter['count'] == total_frames:
            progress = (frame_counter['count'] / total_frames) * 100
            print(f"  フレーム処理中: {frame_counter['count']}/{total_frames} ({progress:.1f}%)", end='\r')

        data = frames_data[frame_idx]
        particles = data['particles']
        
        for idx, circle in enumerate(particle_patches):
            if idx < len(particles):
                p_data = particles[idx]
                circle.center = (p_data['x'], p_data['z'])
                circle.radius = p_data['r']
                circle.set_visible(True)

                x_center = p_data['x']
                z_center = p_data['z']
                radius = p_data['r']
                rotation_angle = p_data['rotation_angle']
                x_end = x_center + radius * np.cos(rotation_angle)
                z_end = z_center + radius * np.sin(rotation_angle)

                rotation_lines[idx].set_data([x_center, x_end], [z_center, z_end])
                rotation_lines[idx].set_visible(True)
            else:
                circle.set_visible(False)
                rotation_lines[idx].set_visible(False)

        time_text_obj.set_text(f"Time: {data['time']:.6f} s")
        return artists_to_update

    def init_frame():
        return artists_to_update

    print("\nアニメーション作成を開始します...")
    print(f"  使用フレーム数: {len(frames_data)}")
    print(f"  出力ファイル: {output_filename}")

    ani = animation.FuncAnimation(
        fig,
        update_frame,
        frames=len(frames_data),
        init_func=init_frame,
        blit=True,
        interval=150,
    )

    try:
        print("\nアニメーション保存中...")
        output_path = Path(output_filename)
        suffix = output_path.suffix.lower()
        final_output = str(output_path)
        ffmpeg_path = shutil.which("ffmpeg")
        
        if suffix == '.gif':
            writer = animation.PillowWriter(fps=fps)
        elif ffmpeg_path:
            matplotlib.rcParams['animation.ffmpeg_path'] = ffmpeg_path
            writer = animation.FFMpegWriter(fps=fps, codec='libx264', extra_args=['-pix_fmt', 'yuv420p'])
        else:
            print("FFmpeg が見つかりません。GIF 出力へフォールバックします。")
            writer = animation.PillowWriter(fps=fps)
            final_output = str(output_path.with_suffix('.gif'))

        ani.save(final_output, writer=writer, dpi=80)
        print(f"\n\nアニメーション保存完了: {final_output}")
        
    except Exception as e:
        print(f"Error during animation creation: {e}")

def parse_arguments():
    parser = argparse.ArgumentParser(description="PEM アニメーション作成ツール (CSV対応版)")
    parser.add_argument("data_file", nargs="?", default=str(DATA_FILE_DEFAULT), help="解析するデータファイル (particles.csv)")
    parser.add_argument("output_file", nargs="?", default="pem_animation.mp4", help="出力ファイル名")
    parser.add_argument("walls_file", nargs="?", default=str(WALLS_FILE_DEFAULT), help="斜面壁データファイル")
    parser.add_argument("--frame-step", type=int, default=1, help="読み込みステップ間隔")
    parser.add_argument("--max-frames", type=int, default=200, help="最大フレーム数")
    return parser.parse_args()


def main():
    args = parse_arguments()

    data_file = Path(args.data_file).expanduser()
    output_file = Path(args.output_file).expanduser()
    walls_file = Path(args.walls_file).expanduser()

    frame_step = max(1, args.frame_step)
    max_frames = args.max_frames if args.max_frames and args.max_frames > 0 else None

    print("=" * 60)
    print("PEM アニメーション作成ツール (CSV版)")
    print("=" * 60)
    print(f"データファイル: {data_file}")
    print(f"出力ファイル: {output_file}")
    
    walls = read_walls_data(walls_file)
    all_frames_data = read_simulation_data(
        data_file,
        frame_step=frame_step,
        max_frames=max_frames,
    )

    if all_frames_data:
        animate(all_frames_data, str(output_file), walls_data=walls)
        print("\n処理が完了しました!")
    else:
        print(f"\nエラー: {data_file} から有効なデータを読み込めませんでした。")


if __name__ == "__main__":
    main()