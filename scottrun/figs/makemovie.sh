python makeframes.py;
convert -delay 10 -scale 125% pngfiles/*.png crap.mpeg;
open crap.mpeg
