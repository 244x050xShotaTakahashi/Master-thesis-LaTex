# ============================================================================
# 修士論文 LaTeX Makefile (LuaLaTeX専用)
# ============================================================================

# ターゲットファイル（拡張子なし）
TARGET = charged_dust_dem_thesis
LATEX = lualatex

# デフォルトターゲット
all: $(TARGET).pdf

# PDF生成 (LuaLaTeXは直接PDFを出力するためdvipdfmxは不要)
# 相互参照解決のため複数回実行
$(TARGET).pdf: $(TARGET).tex
	$(LATEX) $(TARGET)
	$(LATEX) $(TARGET)
	$(LATEX) $(TARGET)

# latexmkを使用する場合（推奨・自動で依存関係を解決）
latexmk:
	latexmk -lualatex $(TARGET)

# 監視モード
watch:
	latexmk -pvc -lualatex $(TARGET)

# クリーン
clean:
	rm -f $(TARGET).aux $(TARGET).log $(TARGET).toc
	rm -f $(TARGET).lof $(TARGET).lot
	rm -f $(TARGET).out $(TARGET).dvi
	rm -f $(TARGET).bbl $(TARGET).blg
	rm -f $(TARGET).fdb_latexmk $(TARGET).fls
	rm -f $(TARGET).synctex.gz

distclean: clean
	rm -f $(TARGET).pdf

.PHONY: all clean distclean latexmk watch