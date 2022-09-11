\rm -f pdffiles/*.pdf
python makeframes.py;
convert -delay 10 -scale 125% pdffiles/*.pdf crap.mpeg;
open crap.mpeg
