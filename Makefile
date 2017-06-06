clean: 
	rm -r results/*;
	rm -r sources/*;
	rm -r shell_calc.sh.o*;
sync:
	rsync -r ~/NCSM/calc/mcalc/results/vce sources/;
