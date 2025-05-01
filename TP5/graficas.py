import numpy as np
import scipy.io
import scipy.signal as signal
import matplotlib.pyplot as plt
import scienceplots
from control import TransferFunction

datos = scipy.io.loadmat("TP5_caso_d.mat")

sys=TransferFunction(datos["sys"][0,0][0][0][0][0][0][0][0],datos["sys"][0,0][1][0][0][0][0][0][0])
sys_Norm=TransferFunction(datos["sys_Norm"][0,0][0][0][0][0][0][0][0],datos["sys_Norm"][0,0][1][0][0][0][0][0][0])
sysc=TransferFunction(datos["sysc_Norm"][0,0][0][0][0][0][0][0][0],datos["sysc_Norm"][0,0][1][0][0][0][0][0][0])

print("La función de transferencia dada es: ")
print(sys._repr_latex_())

print("La función de transferencia normalizada del sistema es:")
print(sys_Norm._repr_latex_())

print("La función de transferencia identificada normalizada es: ")
print(sysc._repr_latex_())

# Recolección de los datos de la simulación
z=datos.get("z")
z0=datos.get("zo")
ys=datos.get("ys")
offset=datos.get("off_set")
med=datos.get("Med")

t_D=datos.get("t_D")
ts=datos.get("ts")
y_D=datos.get("y_D")
t_sal=datos.get("t_sal")
y_sal=datos.get("y_sal")

y_sal2=datos.get("y_sal2")
t_sal2=datos.get("t_sal2")
y_D2=datos.get("y_D2")
t_D2=datos.get("t_D2")

mag_sys=datos.get("mag_sys")
fase_sys=datos.get("fase_sys")
W=datos.get("W")
mag_sys_id=datos.get("mag_sysid")
fase_sys_id=datos.get("fase_sysid")
W_id=datos.get("Wid")

# Gráfica de la simulación

with plt.style.context("ieee"):
    fig1, axs= plt.subplots(2,1,figsize=(10,5),dpi=100)

fig1.subplots_adjust(right=0.8)
  
axs[0].grid(True)
axs[0].plot(z,".")
axs[0].plot(z0,"+r")
axs[0].plot(ys[0,1+int(offset[0,0]):int(offset[0,0]+med[0,0])],"-b")
axs[0].legend(["Estimada","Mediciones","Real"],loc="center left", bbox_to_anchor=(1.05, 0.5))
axs[0].set_title(r"Ajuste con el orden $b_s/a_s$")
axs[0].set_xlabel("Muestras")

axs[1].grid(True)
axs[1].plot(t_D*ts,y_D,"-b")
axs[1].plot(t_sal*ts,y_sal,"--k")
axs[1].legend(["Real","Identificada"],loc="center left", bbox_to_anchor=(1.05, 0.5))
axs[1].set_title("Desempeño del modelo ajustado")
axs[1].set_xlabel("Tiempo. [s]")
plt.tight_layout()


with plt.style.context("ieee"):
    fig2, axs2= plt.subplots(figsize=(10,5),dpi=100)

fig2.subplots_adjust(right=0.8)
axs2.grid(True)
axs2.plot(t_D2,y_D2,"-b")
axs2.plot(t_sal2,y_sal2,".k")
axs2.set_title("Respuesta a la entrada escalón")
axs2.set_xlabel("Tiempo. [s]")
axs2.legend(["Real","Identificada"],loc="center left",bbox_to_anchor=(1.05, 0.5))
plt.tight_layout()

with plt.style.context("ieee"):
    fig3, axs3= plt.subplots(2,1,figsize=(10,5),dpi=100)

fig3.subplots_adjust(right=0.8)
axs3[0].semilogx(W[0,:],mag_sys[:,0])
axs3[0].semilogx(W_id[0,:],mag_sys_id[:,0])
axs3[0].grid(True)
axs3[0].set_title("Magnitud")
axs3[0].legend(["Real","Identificada"],loc="center left",bbox_to_anchor=(1.05, 0.5))
axs3[0].set_xlabel("Frecuencia [w]")

axs3[1].semilogx(W[0,:],fase_sys[:,0])
axs3[1].semilogx(W_id[0,:],fase_sys_id[:,0])
axs3[1].grid(True)
axs3[1].set_title("Fase")
axs3[1].legend(["Real","Identificada"],loc="center left",bbox_to_anchor=(1.05, 0.5))
axs3[1].set_xlabel("Frecuencia [w]")

plt.tight_layout()
plt.show()
