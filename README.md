# 修士論文 LaTeX 環境

帯電ダストダイナミクス解析のための個別要素法（DEM）シミュレーションモデルの開発に関する修士論文を執筆するためのLaTeX環境です。

## プロジェクト概要

このプロジェクトは、帯電した粉体の挙動解析に関する研究の成果をまとめる修士論文の執筆環境です。主に以下の2つのコンポーネントで構成されています：

- **thesis/**: 修士論文のLaTeXソースコードとコンパイル環境
- **DEM-valid/**: DEMシミュレーションプログラム（[GitHub リポジトリ](https://github.com/244x050xShotaTakahashi/DEM-valid)からクローン）

## ディレクトリ構成

```
DEM-valid-local/
├── thesis/                    # 修士論文のLaTeXソースコード
│   ├── charged_dust_dem_thesis.tex  # メインのLaTeXファイル
│   ├── Makefile              # コンパイル用Makefile
│   ├── figures/              # 図を格納するディレクトリ
│   └── README.md             # 詳細な使い方（thesis/README.mdを参照）
│
├── DEM-valid/                # DEMシミュレーションプログラム
│   ├── src/                  # Fortranソースコード
│   ├── scripts/              # 解析・可視化スクリプト
│   ├── inputs/               # シミュレーション入力ファイル
│   └── README.md             # DEMプログラムの詳細
│
└── README.md                 # このファイル
```

## 論文のコンパイル方法

論文をPDFにコンパイルするには、`thesis/`ディレクトリで以下のいずれかのコマンドを実行してください：

### 推奨方法: Makefile を使用

```bash
cd thesis
make                # LuaLaTeXでPDFを生成
make lualatex       # LuaLaTeXを明示的に使用
make clean          # 中間ファイルを削除
```

### latexmk を使用

```bash
cd thesis
latexmk -lualatex charged_dust_dem_thesis    # 通常コンパイル
latexmk -pvc -lualatex charged_dust_dem_thesis  # 監視モード（自動コンパイル）
```

## 必要な環境

- **LaTeX環境**: TeX Live または MacTeX
  - macOSの場合: `brew install --cask mactex`
- **必須パッケージ**: `luatexja`, `amsmath`, `graphicx`, `hyperref`, `siunitx`, `physics` など
- **推奨エディタ**: VS Code + LaTeX Workshop 拡張機能、TeXShop、TeXstudio など

## 詳細情報

論文のコンパイル方法やLaTeXの使い方の詳細については、[thesis/README.md](thesis/README.md)を参照してください。

DEMシミュレーションプログラムの使用方法については、[DEM-valid/README.md](DEM-valid/README.md)を参照してください。

## 論文テーマ

**タイトル**: 帯電ダストダイナミクス解析のための個別要素法シミュレーションモデルの開発

本研究では、静電気を帯びた粉体の挙動解析のためのDEMプログラムを開発し、その妥当性を検証しています。特に、帯電分布やクーロン力のカットオフ距離が粒子堆積物の安息角に与える影響を評価しています。
