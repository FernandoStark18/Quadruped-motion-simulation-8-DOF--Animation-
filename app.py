# Espericueta Fonseca Fernando Simón
# Mecatrónica 6to 6
# Cinemática de robots Proyecto final
# Análisis cinemático por la convención Denavit-Hartenberg de un cuadrúpedo 8 GDL
# Animación

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

fig, ax = plt.subplots()
ax = plt.axes(projection = "3d")
#Función para rotar en el eje x
def matriz_rotacion_x(grados):
	rad = grados/180*np.pi
	#Matriz de rotación
	rotacion=np.array([[1,0,0,0],
					   [0,np.cos(rad),-np.sin(rad),0],
					   [0,np.sin(rad),np.cos(rad),0],
					   [0,0,0,1]])
	return rotacion
#Función para rotar en el eje z
def matriz_rotacion_z(grados):
	rad = grados/180*np.pi
	#Matriz de rotacuón
	rotacion=np.array([[np.cos(rad),-np.sin(rad),0,0],
					   [np.sin(rad),np.cos(rad),0,0],
					   [0,0,1,0],
					   [0,0,0,1]])
	return rotacion
#Función para trasladar en el eje x
def matriz_traslacion_x(x):
	traslacion = np.array([[1,0,0,x],
						   [0,1,0,0],
						   [0,0,1,0],	   
						   [0,0,0,1]])
	return traslacion
#Función para trasladar en el eje z
def matriz_traslacion_z(z):
	traslacion = np.array([[1,0,0,0],
						   [0,1,0,0],
						   [0,0,1,z],	   
						   [0,0,0,1]])
	return traslacion
# Función para la configuración gráfica
def configuracion_grafica():
	# Título del gráfico
	plt.title("Cuadrúpedo 8 GDL")
	# Límites del gráfico
	ax.set_xlim(-10,19)
	ax.set_ylim(-10,19)
	ax.set_zlim(-10,10)
	# Etiquetas para identificar los ejes
	ax.set_xlabel("x")
	ax.set_ylabel("y")
	ax.set_zlabel("z")
	# configuración de la vista
	ax.view_init(elev=50,azim=90)

	# Aristas del cuadrado a dibujar
	x_cuad = [0, 9, 9, 0]
	y_cuad = [0, 0, 9, 9]
	z_cuad = [0, 0, 0, 0]
	# Creamos los vértices
	vertices = [list(zip(x_cuad,y_cuad,z_cuad))]
	# Creamos una malla dentro del perímetro establecido por los vértices
	poly = Poly3DCollection(vertices, alpha=0.8, color = 'k')
	ax.add_collection3d(poly)
# Función para la operación de Denavit-Hartenberg
def DH(theta_i, di, ai, alpha_i):
	MT = matriz_rotacion_z(theta_i)@matriz_traslacion_z(di)@matriz_traslacion_x(ai)@matriz_rotacion_x(alpha_i)
	return MT
# Composición de las matrices de transformación homogénea
def Cuadrupedo(theta_1, d1, a1, alpha_1, 
			   theta_2, d2, a2, alpha_2,
			   theta_3, d3, a3, alpha_3,
			   theta_4, d4, a4, alpha_4,
			   theta_5, d5, a5, alpha_5,
			   theta_6, d6, a6, alpha_6,
			   theta_7, d7, a7, alpha_7,
			   theta_8, d8, a8, alpha_8):
	
	# Pata 1 Frente derecha
	A0 = np.eye(4)
	_0A1 = DH(theta_1, d1, a1, alpha_1)
	_1A2 = DH(theta_2, d2, a2, alpha_2)
	_0A2 = _0A1@_1A2

	# Pata 2 Frente izquierda
	A02 = np.array([[1,0,0,9],
					[0,1,0,0],
					[0,0,1,0],	   
					[0,0,0,1]])
	_0A12 = DH(theta_3, d3, a3, alpha_3)
	_1A22 = DH(theta_4, d4, a4, alpha_4)
	_0A12 = A02@_0A12
	_0A22 = _0A12@_1A22

	# Pata 3 Atrás izquierda
	A03 = np.array([[1,0,0,9],
					[0,1,0,9],
					[0,0,1,0],	   
					[0,0,0,1]])
	_0A13 = DH(theta_5, d5, a5, alpha_5)
	_1A23 = DH(theta_6, d6, a6, alpha_6)
	_0A13 = A03@_0A13
	_0A23 = _0A13@_1A23

	# Pata 4 Atrás derecha
	A04 = np.array([[1,0,0,0],
					[0,1,0,9],
					[0,0,1,0],	   
					[0,0,0,1]])
	_0A14 = DH(theta_7, d7, a7, alpha_7)
	_1A24 = DH(theta_8, d8, a8, alpha_8)
	_0A14 = A04@_0A14
	_0A24 = _0A14@_1A24

	ax.plot3D([A0[0,3],_0A1[0,3]],[A0[1,3],_0A1[1,3]],[A0[2,3],_0A1[2,3]], color = 'red')
	ax.plot3D([_0A1[0,3],_0A2[0,3]],[_0A1[1,3],_0A2[1,3]],[_0A1[2,3],_0A2[2,3]], color = 'green')
	ax.plot3D([A02[0,3],_0A12[0,3]],[A02[1,3],_0A12[1,3]],[A02[2,3],_0A12[2,3]], color = 'red')
	ax.plot3D([_0A12[0,3],_0A22[0,3]],[_0A12[1,3],_0A22[1,3]],[_0A12[2,3],_0A22[2,3]], color = 'green')
	ax.plot3D([A03[0,3],_0A13[0,3]],[A03[1,3],_0A13[1,3]],[A03[2,3],_0A13[2,3]], color = 'red')
	ax.plot3D([_0A13[0,3],_0A23[0,3]],[_0A13[1,3],_0A23[1,3]],[_0A13[2,3],_0A23[2,3]], color = 'green')
	ax.plot3D([A04[0,3],_0A14[0,3]],[A04[1,3],_0A14[1,3]],[A04[2,3],_0A14[2,3]], color = 'red')
	ax.plot3D([_0A14[0,3],_0A24[0,3]],[_0A14[1,3],_0A24[1,3]],[_0A14[2,3],_0A24[2,3]], color = 'green')

# Función para animar al cuadrúpedo
def Cuadrupedo_animado(theta_1, d1, a1, alpha_1, 
					   theta_2, d2, a2, alpha_2,
					   theta_3, d3, a3, alpha_3,
					   theta_4, d4, a4, alpha_4,
					   theta_5, d5, a5, alpha_5,
					   theta_6, d6, a6, alpha_6,
					   theta_7, d7, a7, alpha_7,
					   theta_8, d8, a8, alpha_8):
	
	for i in range(theta_1,60,-5):
		ax.cla()
		configuracion_grafica()
		Cuadrupedo(i+135,0,4,-90,i,0,5,0,
				 theta_3+225,0,4,-90,theta_4,0,5,0,
				 theta_5+315,0,4,-90,theta_6,0,5,0,
				 theta_7+45 ,0,4,-90,theta_8,0,5,0)
		plt.draw()
		plt.pause(1e-3)

	for j in range(theta_5,120,5):
		theta_6-=5
		ax.cla()
		configuracion_grafica()
		Cuadrupedo(i+135,0,4,-90,i,0,5,0,
				 theta_3+225,0,4,-90,theta_4,0,5,0,
				 j+315,0,4,-90,theta_6,0,5,0,
				 theta_7+45 ,0,4,-90,theta_8,0,5,0)
		plt.draw()
		plt.pause(1e-3)

	for k in range(theta_3,60,-5):
		ax.cla()
		configuracion_grafica()
		Cuadrupedo(i+135,0,4,-90,i,0,5,0,
				 k+225,0,4,-90,k,0,5,0,
				 j+315,0,4,-90,theta_6,0,5,0,
				 theta_7+45 ,0,4,-90,theta_8,0,5,0)
		plt.draw()
		plt.pause(1e-3)

	for l in range(theta_7,120,5):
		theta_8 -=5
		ax.cla()
		configuracion_grafica()
		Cuadrupedo(i+135,0,4,-90,i,0,5,0,
				 k+225,0,4,-90,k,0,5,0,
				 j+315,0,4,-90,theta_6,0,5,0,
				 l+45 ,0,4,-90,theta_8,0,5,0)
		plt.draw()
		plt.pause(1e-3)
		theta_5 = 120
	for m in range(theta_1,120,5):
		theta_5 -= 10
		ax.cla()
		configuracion_grafica()
		Cuadrupedo(m+135,0,4,-90,i,0,5,0,
				 k+225,0,4,-90,k,0,5,0,
				 theta_5+315,0,4,-90,theta_6,0,5,0,
				 l+45 ,0,4,-90,theta_8,0,5,0)
		plt.draw()
		plt.pause(1e-3)
		theta_7 = 120
	for n in range(theta_3,120,5):
		theta_7 -= 10
		ax.cla()
		configuracion_grafica()
		Cuadrupedo(m+135,0,4,-90,i,0,5,0,
				 n+225,0,4,-90,n,0,5,0,
				 theta_5+315,0,4,-90,theta_6,0,5,0,
				 theta_7+45 ,0,4,-90,n,0,5,0)
		plt.draw()
		plt.pause(1e-3)

	for o in range(theta_2,120,5):
		ax.cla()
		configuracion_grafica()
		Cuadrupedo(m+135,0,4,-90,o,0,5,0,
				 n+225,0,4,-90,n,0,5,0,
				 theta_5+315,0,4,-90,o,0,5,0,
				 theta_7+45 ,0,4,-90,n,0,5,0)
		plt.draw()
		plt.pause(1e-3)

	for p in range(theta_4-5,60,-5):
		ax.cla()
		configuracion_grafica()
		Cuadrupedo(m+135,0,4,-90,o,0,5,0,
				 n+225,0,4,-90,p,0,5,0,
				 theta_5+315,0,4,-90,o,0,5,0,
				 theta_7+45 ,0,4,-90,p,0,5,0)
		plt.draw()
		plt.pause(1e-3)

	for q in range(theta_1,60,-5):
		theta_5+=5
		ax.cla()
		configuracion_grafica()
		Cuadrupedo(q+135,0,4,-90,o,0,5,0,
				 n+225,0,4,-90,p,0,5,0,
				 theta_5+315,0,4,-90,o,0,5,0,
				 theta_7+45 ,0,4,-90,p,0,5,0)
		plt.draw()
		plt.pause(1e-3)

	for r in range(theta_3,60,-5):
		theta_7+=5
		ax.cla()
		configuracion_grafica()
		Cuadrupedo(q+135,0,4,-90,o,0,5,0,
				 r+225,0,4,-90,p,0,5,0,
				 theta_5+315,0,4,-90,o,0,5,0,
				 theta_7+45 ,0,4,-90,p,0,5,0)
		plt.draw()
		plt.pause(1e-3)

Cuadrupedo_animado(90,0,4,-90,90,0,5,0,
				 90,0,4,-90,90,0,5,0,
				 90,0,4,-90,90,0,5,0,
				 90,0,4,-90,90,0,5,0)
plt.show()
print("Programado por Fernando Espericueta")