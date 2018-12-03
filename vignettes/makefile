LFILE = paretoconv

#all: knith open 
all: alt

knith: $(LFILE).Rmd
	echo "rmarkdown::render('$(LFILE).Rmd',output_file='$(LFILE).html', runtime='static')" | R --no-save -q

knitr: $(LFILE).Rmd
	echo "rmarkdown::render('$(LFILE).Rmd',rmarkdown::md_document(variant='markdown_github'))" | R --no-save -q

open: $(LFILE).html
	xdg-open $(LFILE).html &

alt: knitr_alt knith_alt clean_alt

knitr_alt: $(LFILE).Rmd
	Rscript -e "rmarkdown::render('$(LFILE).Rmd', rmarkdown::md_document(), clean = FALSE)"

knith_alt: $(LFILE).utf8.md
	Rscript -e "rmarkdown::render('$(LFILE).utf8.md', output_format=rmarkdown::html_document (toc=TRUE, toc_float=TRUE, theme='flatly', number_sections=TRUE))"

clean_alt:
	mv $(LFILE).utf8.html $(LFILE).html; \
	rm $(LFILE).md $(LFILE).utf8.md $(LFILE).knit.md

open_alt: $(LFILE).utf8.html
	xdg-open $(LFILE).utf8.html &

clean:
	rm -rf *.html *.png README_cache 
