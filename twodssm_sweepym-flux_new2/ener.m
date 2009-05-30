pl	1	#
		da $1
		erase
		read {t 1 ek 2 eth 3 eg 4 ekx 5 eky 6 px 7 py 8 ag 9 ah 10 vxt 11 st 12 vyt 13 et 14}
		set lvxt=ln(abs(vxt))
		set lst=ln(abs(st))
		set lvyt=ln(abs(vyt))
		set let=ln(abs(et))
		set etot=ek+eth+eg
		#
		window 1 2 1 2
		limits t (ek concat eth concat eg concat etot)
		box
		#connect t etot
		ctype yellow
		connect t ek
		ctype green
		connect t eth
		ctype red
		connect t eg
		ctype white
		xla t \Omega
		yla E
		#
		window 1 2 1 1
		limits t (ag concat ah)
		box
		ltype 0
		connect t ag
		ltype 3
		connect t ah
		ltype 0
		xla t \Omega
		yla \alpha_{eff}
		#
		#apl 1 2 ek "kinetic energy"
		#apl 2 2 eb "magnetic energy"
		#apl 1 1 eg "gravitational energy"
		#apl 2 1 de "energy error"
		#
		window 1 1 1 1
		#
		set atot = ag+ah
		stats atot me si ku
		echo "ave alpha: " $me
		#
		set g = 2.
		set Q = sqrt(g*(g-1.)*eth)/pi
		stats Q me si ku
		echo "ave Q: " $me
		limits t lvxt
		erase
		box
		connect t lvxt
		#
apl	4	#
		window 2 2 $1 $2
		limits t $3
		box
		connect t $3
		xla t
		yla $4
bpl	0	#
		da ener.out
		erase
		read {t 1 ek 2 eth 3 eg 4 de 5}
		set etot=ek+eth+eg
		#
		limits t (ek concat eth)
		ctype white
		box
		ctype yellow
		connect t ek
		ctype green
		connect t eth
		ctype red
		#
		window 1 1 1 1
lpl		da ener.out
		erase
		read {t 1 ek 2 eth 3 eg 4 de 5}
		set etot=ek+eth+eg
		set lek = lg(ek + 1.e-3)
		set leth = lg(eth + 1.e-3)
		set leg = lg(abs(eg)+1.e-3)
		#
		limits t lew
		ticksize 0 0 -1 0
		box
		connect t lew
		#
		#limits t (lek concat leth concat leg)
		#ctype white
		#box
		#ctype yellow
		#connect t lek
		#ctype green
		#connect t leth
		#ctype red
		#ctype cyan
		#connect t leg
		#ctype white
		#
av	1	#
		da ener.out
		read {t 1 ek 2 eth 3 eg 4 ekx 5 eky 6 px 7 py 8 ag 9 ah 10}
		set etot=ek+eth+eg
		#
		set atot = ag+ah if (t > $1)
		stats atot me si ku
		echo "ave alpha: " $me
		#
		set g = 11./7.
		set Q = sqrt(g*(g-1.)*eth)/pi
		stats Q me si ku
		echo "ave Q: " $me
		#
rd	1	#
		da $1
		read {t 1 ek 2 eth 3 eg 4 ekx 5 eky 6 px 7 py 8 ag 9 ah 10}
		set etot=ek+eth+eg
		set lekx = lg(ekx)
		set leky = lg(eky)
