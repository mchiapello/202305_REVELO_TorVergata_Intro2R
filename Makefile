index: index.qmd
	quarto render index.qmd

M1: materials/M1_intro.qmd
	quarto render materials/M1_intro.qmd

M3: materials/M3_introR.qmd
	quarto render materials/M3_introR.qmd

clean:
	rm docs/index.html
	rm -rf docs/materials/
