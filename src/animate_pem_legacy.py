import argparse
import shutil
from pathlib import Path
from typing import List, Optional, Union

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Circle

DATA_FILE_DEFAULT = Path(__file__).resolve().parent.parent / "data" / "graph11.d"
WALLS_FILE_DEFAULT = Path(__file__).resolve().parent.parent / "inputs" / "walls.dat"

def read_walls_data(filename: Union[str, Path] = WALLS_FILE_DEFAULT):
    """walls.dat を読み込み、斜面壁のリストを返す。"""
    walls = []
    try:
        with open(filename, "r") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('!'):
                    continue
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
                    continue
    except FileNotFoundError:
        pass
    return walls

def read_simulation_data(
    filename: Union[str, Path] = DATA_FILE_DEFAULT,
    frame_step: int = 1,
    max_frames: Optional[int] = None,
):
    """graph11.d (旧形式) を読み込み、各タイムステップのデータを辞書のリストとして返す。"""

    if frame_step < 1:
        raise ValueError("frame_step must be >= 1")
    if max_frames is not None and max_frames < 1:
        raise ValueError("max_frames must be >= 1 when provided")

    def _read_floats(token_buffer, fh, n):
        vals = []
        while len(vals) < n:
            if not token_buffer:
                line = fh.readline()
                if not line:
                    raise EOFError("予期せぬ EOF")
                token_buffer.extend(line.split())
            vals.append(float(token_buffer.pop(0)))
        return vals

    def _skip_floats(token_buffer, fh, n):
        count = 0
        while count < n:
            if not token_buffer:
                line = fh.readline()
                if not line:
                    raise EOFError("予期せぬ EOF")
                token_buffer.extend(line.split())
            token_buffer.pop(0)
            count += 1

    frames_data = []
    print("データファイルを読み込み中(Legacy)...")
    try:
        with open(filename, "r") as fh:
            tokens: list[str] = []
            frame_count = 0
            stored_frames = 0
            while True:
                try:
                    # ヘッダー: num_particles, time, container_width, container_height, rmax
                    num_particles_f, time_val, container_width, container_height, rmax_val = _read_floats(tokens, fh, 5)
                    num_particles = int(num_particles_f)
                except EOFError:
                    break

                if num_particles < 0:
                    continue

                take_frame = (frame_count % frame_step == 0)
                particles_data_to_read = num_particles * 3

                if num_particles > 0:
                    if take_frame:
                        particle_values = _read_floats(tokens, fh, particles_data_to_read)
                    else:
                        _skip_floats(tokens, fh, particles_data_to_read)
                else:
                    particle_values = []

                # 速度データ読み飛ばし
                if num_particles > 0:
                    _skip_floats(tokens, fh, particles_data_to_read)

                # 回転角読み込み
                if num_particles > 0:
                    if take_frame:
                        rotation_angles = _read_floats(tokens, fh, num_particles)
                    else:
                        _skip_floats(tokens, fh, num_particles)
                else:
                    rotation_angles = []

                if take_frame:
                    particles = [
                        {
                            "x": particle_values[i * 3 + 0],
                            "z": particle_values[i * 3 + 1],
                            "r": particle_values[i * 3 + 2],
                            "rotation_angle": rotation_angles[i] if i < len(rotation_angles) else 0.0,
                        }
                        for i in range(num_particles)
                    ]

                    frames_data.append(
                        {
                            "time": time_val,
                            "num_particles": num_particles,
                            "container_width": container_width,
                            "container_height": container_height,
                            "particles": particles,
                        }
                    )
                    stored_frames += 1
                    if max_frames is not None and stored_frames >= max_frames:
                        print(f"\r  {len(frames_data)} フレームの読み込み完了 (制限到達)!        ")
                        break

                frame_count += 1
                if frame_count % 100 == 0:
                    print(f"  {frame_count} フレーム読み込み完了...", end='\r')
    except FileNotFoundError:
        print(f"エラー: ファイル '{filename}' が見つかりません。")
        return None
    except Exception as e:
        print(f"ファイル読み込み中にエラーが発生しました: {e}")
        return None
    
    print(f"\r  {len(frames_data)} フレームの読み込み完了!        ")
    return frames_data

def animate(frames_data, output_filename="pem_animation_legacy.mp4", walls_data=None, fps: int = 10):
    if not frames_data:
        print("アニメーションするデータがありません。")
        return

    import matplotlib
    matplotlib.rcParams['font.family'] = 'DejaVu Sans'

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel("X coordinate")
    ax.set_ylabel("Z coordinate")
    fig.suptitle("PEM Validation (Legacy Format)", fontsize=12)

    if walls_data is None:
        walls_data = []

    all_widths = [frame['container_width'] for frame in frames_data]
    max_container_width = max(all_widths) if all_widths else 1.0
    
    # 高さ方向の最大値
    max_z_overall = 0.0
    for frame_data in frames_data:
        if frame_data['particles']:
            current_max_z = max(p['z'] + p['r'] for p in frame_data['particles'])
            if current_max_z > max_z_overall:
                max_z_overall = current_max_z
    
    # コンテナ高さも含めて表示範囲決定
    container_heights = [frame.get('container_height', 0.0) for frame in frames_data]
    max_container_height = max(container_heights) if container_heights else 0.0
    
    display_height = max(max_z_overall, max_container_height)
    if display_height <= 0: display_height = 1.0
    
    # 少し余裕を持たせる
    display_height *= 1.2

    x_margin = 0.1 * max_container_width
    ax.set_xlim(-x_margin, max_container_width + x_margin)
    ax.set_ylim(-0.1 * display_height, display_height)

    time_text_obj = ax.text(0.05, 0.95, '', transform=ax.transAxes, ha="left", va="top", fontsize=10)

    # 壁
    left_wall = Line2D([], [], color='black', lw=2)
    bottom_wall = Line2D([], [], color='black', lw=2)
    right_wall = Line2D([], [], color='black', lw=2)
    top_wall = Line2D([], [], color='black', lw=2)
    for wall_line in (left_wall, bottom_wall, right_wall, top_wall):
        ax.add_line(wall_line)
    
    if walls_data:
        for wall in walls_data:
            line = Line2D(
                [wall['x_start'], wall['x_end']],
                [wall['z_start'], wall['z_end']],
                color='black', lw=2.5, alpha=0.8
            )
            ax.add_line(line)

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

    artists_to_update = [time_text_obj, left_wall, bottom_wall, right_wall, top_wall]
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
        cw = data['container_width']
        ch = data.get('container_height', 0.0)
        
        wall_h = ch if ch > 0.0 else display_height
        left_wall.set_data([0, 0], [0, wall_h])
        bottom_wall.set_data([0, cw], [0, 0])
        right_wall.set_data([cw, cw], [0, wall_h])
        
        if ch > 0.0:
            top_wall.set_data([0, cw], [ch, ch])
            top_wall.set_visible(True)
        else:
            top_wall.set_visible(False)

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

    ani = animation.FuncAnimation(
        fig, update_frame, frames=len(frames_data), init_func=init_frame, blit=True, interval=150
    )

    print(f"\nアニメーション保存中: {output_filename}")
    output_path = Path(output_filename)
    suffix = output_path.suffix.lower()
    
    ffmpeg_path = shutil.which("ffmpeg")
    if suffix == '.gif':
        writer = animation.PillowWriter(fps=fps)
    elif ffmpeg_path:
        matplotlib.rcParams['animation.ffmpeg_path'] = ffmpeg_path
        writer = animation.FFMpegWriter(fps=fps, codec='libx264', extra_args=['-pix_fmt', 'yuv420p'])
    else:
        writer = animation.PillowWriter(fps=fps)
        output_path = output_path.with_suffix('.gif')
        
    try:
        ani.save(str(output_path), writer=writer, dpi=80)
        print(f"保存完了: {output_path}")
    except Exception as e:
        print(f"保存エラー: {e}")

def main():
    parser = argparse.ArgumentParser(description="PEM アニメーション作成ツール (Legacy graph11.d対応版)")
    parser.add_argument("data_file", nargs="?", default=str(DATA_FILE_DEFAULT), help="解析するデータファイル (graph11.d)")
    parser.add_argument("output_file", nargs="?", default="pem_animation_legacy.mp4", help="出力ファイル名")
    parser.add_argument("walls_file", nargs="?", default=str(WALLS_FILE_DEFAULT), help="斜面壁データファイル")
    parser.add_argument("--frame-step", type=int, default=1, help="読み込みステップ間隔")
    parser.add_argument("--max-frames", type=int, default=200, help="最大フレーム数")
    args = parser.parse_args()

    data_file = Path(args.data_file).expanduser()
    output_file = Path(args.output_file).expanduser()
    walls_file = Path(args.walls_file).expanduser()

    print("=" * 60)
    print("PEM アニメーション作成 (Legacy)")
    print(f"Input: {data_file}")
    print(f"Output: {output_file}")
    print("=" * 60)

    walls = read_walls_data(walls_file)
    all_frames_data = read_simulation_data(
        data_file,
        frame_step=max(1, args.frame_step),
        max_frames=args.max_frames if args.max_frames and args.max_frames > 0 else None,
    )

    if all_frames_data:
        animate(all_frames_data, str(output_file), walls_data=walls)
    else:
        print("データ読み込み失敗")

if __name__ == "__main__":
    main()