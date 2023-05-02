index: index.qmd
	quarto render index.qmd

M1: materials/M1_intro.qmd
	quarto render materials/M1_intro.qmd

clean:
	rm docs/index.html
	rm -rf docs/materials/
