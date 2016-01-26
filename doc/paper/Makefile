TARGET=paper
all:
	pdflatex $(TARGET)
	bibtex $(TARGET).aux
	pdflatex $(TARGET)

clean:
	rm -f $(TARGET).bbl  $(TARGET).log  $(TARGET).pdf *.aux  $(TARGET).blg  $(TARGET).out
