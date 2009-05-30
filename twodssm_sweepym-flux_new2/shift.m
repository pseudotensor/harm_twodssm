pl  2	#
	data $1
	lines 1 64
	read {j 1 si1 2 vi1 3 so1 4 vo1 5}
	lines 65 128
	read {sti 2 si2 3 vi2 4 so2 5 vo2 6}
	limits j $2o1
	erase
	ctype white
	box
	connect j $2o1
	ctype yellow
	connect j $2i2
	ctype red
	limits j sti
	connect j sti
