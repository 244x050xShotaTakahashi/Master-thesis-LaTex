# 修士論文 LaTeX テンプレート

帯電ダストダイナミクス解析のための個別要素法シミュレーションモデルの開発に関する修士論文のLaTeXテンプレートです。

## ファイル構成

```
thesis/
├── charged_dust_dem_thesis.tex  # メインのLaTeXファイル
├── Makefile                      # コンパイル用Makefile
├── README.md                     # このファイル
├── figures/                      # 図を格納するディレクトリ
│   └── .gitkeep
└── latexmkrc                     # latexmk設定ファイル
```

## コンパイル方法

### 方法1: latexmk（推奨）

```bash
# PDFを生成
latexmk -pdf charged_dust_dem_thesis

# LuaLaTeXを使用（日本語対応が良好）
latexmk -lualatex charged_dust_dem_thesis

# 監視モード（ファイル変更時に自動コンパイル）
latexmk -pvc -pdf charged_dust_dem_thesis
```

### 方法2: Makefile

```bash

# LuaLaTeXを使用
make 

# あるいは
make lualatex

# クリーン（中間ファイル削除）
make clean

# 完全クリーン（PDFも削除）
make distclean
```

### 方法3: 手動コンパイル

```bash
# LuaLaTeX（推奨）
lualatex charged_dust_dem_thesis
lualatex charged_dust_dem_thesis
lualatex charged_dust_dem_thesis

```

## 必要なパッケージ

このテンプレートを使用するには、以下のLaTeXパッケージが必要です：

- ltjsreport（日本語レポートクラス）
- amsmath, amssymb, amsthm（数式）
- graphicx（図の挿入）
- booktabs, tabularx, multirow（表）
- listings（ソースコード）
- algorithm, algorithmic（アルゴリズム）
- hyperref（ハイパーリンク）
- siunitx（SI単位系）
- physics（物理記号）
- subcaption（サブキャプション）
- titlesec（章の余白調整）
- geometry（ページレイアウト）

### TeX Liveのインストール

macOSの場合：
```bash
# Homebrewを使用
brew install --cask mactex

# または BasicTeX（軽量版）
brew install --cask basictex
```

## 図の追加方法

1. `figures/` ディレクトリに図ファイル（PDF, PNG, JPG等）を配置
2. 本文中で以下のように参照：

```latex
\begin{figure}[htbp]
\centering
\includegraphics[width=0.8\textwidth]{figures/figure_name.pdf}
\caption{図のキャプション}
\label{fig:label}
\end{figure}
```

## 編集時の注意点

1. **章の追加**: `\chapter{章タイトル}` で新しい章を追加
2. **節の追加**: `\section{節タイトル}` で新しい節を追加
3. **数式**: `\begin{equation}...\end{equation}` で番号付き数式
4. **参考文献**: `\bibitem{key}` で参考文献を追加、`\cite{key}` で引用
5. **相互参照**: `\label{xxx}` でラベルを付け、`\ref{xxx}` で参照

## カスタマイズ

### タイトル・著者情報の変更

`charged_dust_dem_thesis.tex` の以下の部分を編集してください：

```latex
\title{...}
\author{...}
\date{...}
```

### ページレイアウトの変更

`charged_dust_dem_thesis.tex` の `geometry` パッケージの設定を変更：

```latex
\usepackage[top=30mm, bottom=30mm, inner=30mm, outer=25mm]{geometry}
```

## トラブルシューティング

### 日本語が表示されない

- LuaLaTeXを使用してください
- または、platex + dvipdfmx の組み合わせを使用

### 図が表示されない

- 図ファイルのパスが正しいか確認
- PDF形式の図を推奨（ベクタ形式で高品質）

### コンパイルエラー

- 中間ファイルを削除して再コンパイル：`make clean && make`
- エラーメッセージを確認してパッケージの不足を確認











