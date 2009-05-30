#
for resx in 64 128 256 512 1024
do
	let resy=$((2*$resx))
	echo $resy
	sed "s/RESX/$resx/g;s/RESY/$resy/g" < decs.base > decs.h
	make twodssm
	twodssm
	mv ener.out ener.out.$resx
	mkdir run.$resx
	mv images/im0* run.$resx
	mv dump0* run.$resx
done
