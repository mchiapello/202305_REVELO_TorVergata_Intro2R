index: index.qmd
	quarto render index.qmd

M1: materials/M1_intro.qmd
	quarto render materials/M1_intro.qmd
	open docs/materials/M1_intro.html -a /Applications/Safari.app/

clean:
	rm docs/index.html
	rm -rf docs/materials/
