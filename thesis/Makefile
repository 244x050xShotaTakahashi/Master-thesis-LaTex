# ============================================================================
# 修士論文 LaTeX Makefile
# ============================================================================

# ターゲットファイル
TARGET = charged_dust_dem_thesis
LATEX = lualatex
BIBTEX = pbibtex
DVIPDF = dvipdfmx

# デフォルトターゲット
all: $(TARGET).pdf

# PDF生成 (LuaLaTeXは直接PDFを出力するためdvipdfmxは不要)
$(TARGET).pdf: $(TARGET).tex
	$(LATEX) $(TARGET)
	$(LATEX) $(TARGET)

# 旧設定（コメントアウト）
#$(TARGET).dvi: $(TARGET).tex
#	$(LATEX) $(TARGET)
#	$(LATEX) $(TARGET)
#	$(LATEX) $(TARGET)
#
#$(TARGET).pdf: $(TARGET).dvi
#	$(DVIPDF) $(TARGET)

# LuaLaTeXを使用する場合（推奨）
lualatex:
	lualatex $(TARGET)
	lualatex $(TARGET)
	lualatex $(TARGET)

# pdfLaTeXを使用する場合
pdflatex:
	pdflatex $(TARGET)
	pdflatex $(TARGET)
	pdflatex $(TARGET)

# latexmkを使用する場合（自動で依存関係を解決）
latexmk:
	latexmk -pdf $(TARGET)

# latexmk + LuaLaTeX
latexmk-lua:
	latexmk -lualatex $(TARGET)

# クリーン（中間ファイルの削除）
clean:
	rm -f $(TARGET).aux $(TARGET).log $(TARGET).toc
	rm -f $(TARGET).lof $(TARGET).lot
	rm -f $(TARGET).out $(TARGET).dvi
	rm -f $(TARGET).bbl $(TARGET).blg
	rm -f $(TARGET).fdb_latexmk $(TARGET).fls
	rm -f $(TARGET).synctex.gz

# 完全クリーン（PDFも削除）
distclean: clean
	rm -f $(TARGET).pdf

# 監視モード（ファイル変更時に自動コンパイル）
watch:
	latexmk -pvc -pdf $(TARGET)

.PHONY: all clean distclean lualatex pdflatex latexmk latexmk-lua watch






