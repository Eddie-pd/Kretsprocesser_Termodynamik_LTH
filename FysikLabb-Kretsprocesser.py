#Program för att visualisera uppg i Labb om kretsprocesser Fysik-Termodynamik.
import numpy as np
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as polynomial
import scipy.optimize as optimize
import pandas


data = pandas.read_csv("Data_Heat_Pump_Del_2.txt", sep='\t', decimal=',', skiprows=2, header=0 )
data.columns = ["Tid", "Count", "T_h", "T_c"]


tiden = np.arange(0, 1360, 10) #Behövde skapa denna då tidskolumnen blev fel i tid arrayen.


tid = data["Tid"].to_numpy() ##Konvertera till numpy array. #riktig count.
count = data["Count"].to_numpy() ##Konvertera till numpy array #Riktig t_h
t_h = data["T_h"].to_numpy()  #Riktig t_c
t_c = data["T_c"].to_numpy() 
w = np.cumsum(tid*4500) ##Addera 4500 för varje 1:a

#print(data)
t_hny = np.array(count+273.15) #Omvandla till Kelvin
#e_c = np.transpose(np.array([t_hny/(count-t_h)])) #Beräkna Carnot-värmefaktor
e_c = np.transpose([t_hny/(count-t_h)])
#funktion_e_c = np.poly1d(np.polyfit(tiden, e_c, 5, rcond=None, full=False, w=None, cov=False))
#polyval_e_c = np.polyval(funktion_e_c)
#print(e_c)
t_cny = np.array(t_h + 273.15) #Omvandla till KELvin

funktion_W = np.poly1d(np.polyfit(tiden, w, 1, rcond=None, full=False, w=None, cov=False))
funktion_Th = np.poly1d(np.polyfit(tiden, t_hny*41800, 5, rcond=None, full=False, w=None, cov=False))
derivative_Th = funktion_Th.deriv()
derivative_W = funktion_W.deriv()
#print(derivative_Th, derivative_W)
#lutning_Th = np.polyfit(tiden, t_hny, 5, rcond=None, full=False, w=None, cov=False) #Beräknar lutning på Qh
#lutning_W = np.polyfit(tiden, w, 1, rcond=None, full=False, w=None, cov=False) #Beräknar lutning på W, och nej denna lutning ska inte vara 0, då vi räknar lutning på totala arbete gjort inte lutningen på arbete per sekund.
#Vi ska använda polyval för att kunna plotta w och lutning th, och sen dela dem o få riktigt värmefaktor.
#print(lutning_riktig_Vf)
#Wpolyval = np.polyval(funktion_W, tiden) #Skapar en array med värden som uppfyller en funktion med den lutning och startvärde som står i lutning_W.
#Om man plottar polyvalen ovan så får man en rät-linje version av det man hade fått om man bara plottade alla w värden direkt(en massa wiggly lines), med det är ej slutmålet här.
#Lutningen för den verkliga värmefaktorn är dTh / dw vilket alltså blir lutning_Th[0]/lutning_W[0]. Sen behövs oxå ett startvärde, sen kan det göras till en np array som kan plottas
#THpolyval = np.polyval(funktion_Th, tiden)

TH_derivative_polyval = np.polyval(derivative_Th, tiden)
W_derivative_polyval = np.polyval(derivative_W, tiden)
#dTH_polyval = np.polyval(np.array([[lutning_Th[0]*2], [lutning_Th[1]]]), tiden)
Vf_riktig = np.asarray(TH_derivative_polyval/W_derivative_polyval)

#Tänker försöka skapa en array a som har samma format som e_c där alla värdena från Vf_riktig ska finnas.
a = [0]*len(e_c)

for i in range(len(Vf_riktig)):
   a[i] =  e_c[i][0]
  
#print(a)
Vf_riktig.reshape((136, 1))

#Vf_kvot = np.transpose(np.transpose(Vf_riktig)/e_c)
#Vf_kvot = np.divide(Vf_riktig, e_c)
Vf_kvot = Vf_riktig/a 
#Vf_kvot = np.array(np.transpose(Vf_kvot))
#print(W_derivative_polyval)
#print(TH_derivative_polyval)
#print("Carnot: ", e_c)
print("Vf kvot: " , Vf_kvot)
print("Vf_riktig", Vf_riktig)
#print(Vf_kvot)

plt.plot(tiden, Vf_kvot, "hotpink")
plt.title(r"Kvoten mellan verklig och ideal värmefaktor över tid")
plt.xlabel("Tid [s]")
plt.ylabel("Kvot av värmefaktor")
plt.savefig("Vf_kvot_riktig_graf")
plt.show()


#plt.plot(tiden, e_c, "hotpink")
plt.plot(tiden, Vf_riktig, "green")
plt.title(r"Verklig Värmefaktor över tid")
plt.xlabel("Tid [s]")
plt.ylabel("Värmefaktor")
plt.savefig("Vf_riktig_graf")
plt.show()

plt.plot(tiden, e_c, "hotpink")
plt.title(r"Carnots värmefaktor $Vf_C$ per tidsenhet $t$")
plt.xlabel("Tid [s]")
plt.ylabel("Värmefaktor (Carnots)")
plt.savefig("ec_graf")
plt.show()


plt.plot(tiden, t_hny, "hotpink")
plt.title(r"Temperatur varm reservoar per tidsenhet")
plt.xlabel("Tid [s]")
plt.ylabel("Th[K]")
plt.savefig("Th_graf")
plt.show()

plt.plot(tiden, t_cny, "hotpink")
plt.title(r"Temperatur kall reservoar per tidsenhet")
plt.xlabel("Tid [s]")
plt.ylabel("Tc/Tid")
plt.savefig("Tc_graf")
plt.show()

plt.plot(tiden, w, "hotpink")
plt.title(r"Arbete in i kompressor per tidsenhet")
plt.xlabel("Tid [s]")
plt.ylabel("W[J]")
plt.savefig("W_graf")
plt.show()


