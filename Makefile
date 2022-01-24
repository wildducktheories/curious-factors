publish:
	jupyter nbconvert README.ipynb --LatexPreprocessor.enabled True --to html --stdout > README.html

