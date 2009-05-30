#
for res in 32 64 128 256 512 1024
do
	sed "s/RES/$res/g" < decs.base > decs.h
	make twodssm
	twodssm
	mv ener.out ener.out.$res
	mkdir run.$res
	mv images/im0* run.$res
	mv dump0* run.$res
done
