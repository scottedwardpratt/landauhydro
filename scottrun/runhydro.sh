make hydro; rm -r output/*.dat; hydro; cd figs; python xslice_animate.py; open -a Safari density.png; cd ~/git/landauhydro/scottrun