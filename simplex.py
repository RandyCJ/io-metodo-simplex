import sys 
import os
import sympy
from iteration_utilities import duplicates
from sympy import *
from sympy.core.rules import Transform

class Matriz:

    """ Clase matriz, contiene todos los métodos necesarios para realizar las operaciones del método simplex,
       y contiene también los atributos necesarios del método simplex
    """

    matriz = [] # La matriz con la fila 0 con las variables y con la columna 0 con las variables básicas
    variables_basicas = [] # Contiene las variables básicas en la columna 0 en un determinado punto
    columna_pivote = () # Contiene la posición de la variable básica entrante y la misma variable entrante de una iteración
    fila_pivote = () # Contiene la posición de la variable básica saliente y la misma variable saliente de una iteración
    indice_pivote = ()
    pivote = 1
    acotada = False
    degenerada = False
    var_degeneradas = [] # Contiene las variables que son degeneradas si existen
    soluciones_multiples = False
    nom_archivo = ""
    U = 0
    FEV = []
    dual = False
    CONST_M = 0
    var_artificiales = []
    is_max = True

    def __init__(self, dict_datos) -> None:
        self.definir_ecuaciones(dict_datos)

    def crear_matriz(self, matriz, variables_basicas, cant_variables):
        """ Crea la matriz junto con el encabezado de las variables y los números y valores de cada una
            E: una matriz con solo los valores, las variables que van en la columna 0 y la cantidad de variables que deben haber
            S: N/A
        """
        encabezado = ["VB"]

        for i in range(cant_variables):
            encabezado.append("X" + str(i+1))
        encabezado.append("LD")
        
        i = 0
        nueva_matriz = []

        for fila in matriz:
            nueva_fila = []
            if i == 0:
                if self.is_max:
                    nueva_fila.append("U")
                else:
                    nueva_fila.append("-U")
                nueva_fila += fila
                nueva_matriz.append(nueva_fila)
                i += 1
            else:
                nueva_fila.append("X" + str(variables_basicas[i-1]))
                nueva_fila += fila
                nueva_matriz.append(nueva_fila)
                i += 1

        self.matriz = [encabezado] + nueva_matriz

    def definir_ecuaciones(self, diccionario_datos):
        """ Convierte el diccionario de datos a las ecuaciones que se utilizaran para las tablas
            E: Recibe el diccionario de datos con los datos que se recolectaron al leer el archivo
            S: N/A
        """
        if diccionario_datos["metodo"] == 1:     #granm
            self.definir_ecuaciones_granm(diccionario_datos)
            return
        elif diccionario_datos["metodo"] == 3:   #dual
            self.dual = True
            diccionario_datos = self.acomodar_diccionario(diccionario_datos)
            print(diccionario_datos)

        diccionario_datos["fun_ob"] = [-x for x in diccionario_datos["fun_ob"]] #Se vuelven negativos todos los números
        diccionario_datos["fun_ob"] += [Rational(0)] * (diccionario_datos["num_rest"] + 1) #Agrega los ceros dependiendo de la cantidad de restricciones

        for rest in diccionario_datos["rest"]: #Coloca los ceros y unos en las restricciones
            tmp = rest.pop(-1)
            rest += [Rational(0)] * diccionario_datos["num_rest"]
            rest += [Rational(tmp)]

        i = diccionario_datos["num_var"]
        for rest in diccionario_datos["rest"]:
            rest[i] = 1
            i += 1

        var_basicas = [x for x in range(diccionario_datos["num_var"]+1, len(diccionario_datos["fun_ob"]))]
        self.crear_matriz([diccionario_datos["fun_ob"]] + diccionario_datos["rest"], var_basicas, diccionario_datos["num_var"] + diccionario_datos["num_rest"])

    def encontrar_entrante(self):
        """ Encuentra cual es la columna con la variable entrante, esto viendo cual tiene el menor valor en U
            E: N/A
            S: N/A
        """
        if self.soluciones_multiples: #Si es multiple no hace falta sacar columna entrante
            return 0
        
        fila = self.matriz[1][1:-1]
        for valor in fila:
            if isinstance(valor, sympy.Basic):
                fila[fila.index(valor)] = valor.subs({self.CONST_M:1000})
    
        tmp = sorted(fila)
        idx = fila.index(tmp[0]) + 1
        self.columna_pivote = (self.matriz[0][idx], idx)

    def encontrar_saliente(self):
        """ Encuentra cual es la columna con la variable entrante, esto calculando el pivote
            E: N/A
            S: N/A
        """

        lista_divisiones = []
        lista_degenerados = []
        self.var_degeneradas = []
        self.acotada = True
        vb_saliente = ""

        for i in self.matriz[2:]:
            if i[self.columna_pivote[1]] > 0 and i[-1] >= 0:
                self.acotada = False
                lista_divisiones += [(i[-1]/i[self.columna_pivote[1]], i[0])]
                lista_degenerados += [i[-1]/i[self.columna_pivote[1]]] #Guardamos todos los resultados para al final ver si hay degenerados

        if self.acotada: # Si ninguna division es valida sería no acotada
            msj_acotada = 'La columna ' + str(self.columna_pivote[1]+1)
            msj_acotada += ' con la VB entrante ' + self.columna_pivote[0]
            msj_acotada += ' no tiene números mayores a 0'
            msj_acotada += "\nPor ello la solución es no acotada"
            escribir_archivo(self.nom_archivo,"\n"+msj_acotada)
            print(msj_acotada)
            quit()

        lista_divisiones.sort()
        vb_saliente = lista_divisiones[0][1]
        ind_saliente = 0

        for fila in self.matriz: #Calculamos el índice de la variable entrante
            if fila[0] != vb_saliente:
                ind_saliente += 1
                continue
            break

        self.fila_pivote = (vb_saliente, ind_saliente)
         
        if lista_divisiones[0][0] in list(duplicates(lista_degenerados)): #Revisamos que no hayan divisiones iguales que serían degeneradas
            self.degenerada = True
            i = 0
            self.var_degeneradas += [lista_divisiones[0][0]]
            while i < len(lista_divisiones):
                if lista_divisiones[i][0] == lista_divisiones[0][0]:
                    self.var_degeneradas += [lista_divisiones[i][1]]
                i += 1
        if self.degenerada and self.var_degeneradas != []:
            print(self.degeneradas_a_texto())

    
    def encontrar_basicas(self):
        """ Guarda las variables básicas en la columna 0 de cada fila
            E: N/A
            S: N/A
        """
        var_basicas = []
        for fila in self.matriz[2:]:
            var_basicas.append(fila[0])
        self.variables_basicas = var_basicas

    def encontrar_FEV(self):
        """ Encuentra la FEV de la matriz actual
            E: N/A
            S: N/A
        """
        self.U = self.matriz[1][-1]
        columna = 1
        self.FEV = []
        self.encontrar_basicas()

        for i in self.matriz[1][1:]:
            if i == 0:
                if not(self.matriz[0][columna] in self.variables_basicas) and columna < len(self.matriz[0])-1:
                    self.columna_pivote = (self.matriz[0][columna], columna)
                    self.FEV += [0]
                else:
                    fila = 1
                    while (columna < len(self.matriz[0]) and fila < len(self.matriz)):
                        if self.matriz[fila][columna] == 1:
                            self.FEV += [self.matriz[fila][-1]]   
                        
                        fila += 1

            elif columna < len(self.matriz[0])-1:
                self.FEV += [0]

            columna += 1
    
    def obtener_columna_pivote(self):
        """ Almacena la columna pivote donde está el pivote antes de convertirla en ceros
            E: N/A
            S: N/A
        """
         
        fila = 1
        columna = []
        while (fila < len(self.matriz)):
            columna += [-self.matriz[fila][self.columna_pivote[1]]]
            fila += 1
        self.columna_pivote = (self.columna_pivote[0], self.columna_pivote[1], columna)

    def modificar_linea_pivote(self):
        """ Modifica la línea y columna donde se encuentra el pivote
            La linea la divide entre el pivote, y la columna convierte todo en ceros.
            E: N/A
            S: N/A
        """
        #Primero se convierten todos los números de la columna donde está el pivote a 0
        self.pivote = self.matriz[self.fila_pivote[1]][self.columna_pivote[1]]
        i = 1
        while i < len(self.matriz):
            self.matriz[i][self.columna_pivote[1]] = Rational(0)
            i += 1

        #Después se divide toda la linea entre el pivote
        i = 1
        self.matriz[self.fila_pivote[1]][0] = self.columna_pivote[0] #Sale VB Saliente y entra la VB Entrante
        while i < len(self.matriz[0]):
            self.matriz[self.fila_pivote[1]][i] = self.matriz[self.fila_pivote[1]][i]/self.pivote
            i += 1
        self.matriz[self.fila_pivote[1]][self.columna_pivote[1]] = Rational(1)

    def iteracion(self):
        """ Hace las operaciones en cada fila para convertir la columna pivote en ceros
            E: N/A
            S: N/A
        """
        num_fila = 0
        for fila in self.matriz:
            if num_fila != self.fila_pivote[1] and num_fila != 0 :
                indice = 0
                for num in fila[1:]:
                    if indice != self.columna_pivote[1] - 1 and indice <= len(self.matriz[0]) + 1:
                        self.matriz[num_fila][indice + 1] = num + self.columna_pivote[2][num_fila - 1] * self.matriz[self.fila_pivote[1]][indice + 1]
                    indice += 1                
            num_fila += 1
    
    def iterar(self):
        """ Función encargada de llamar a las operaciones necesarias para realizar una iteración
            E: N/A
            S: N/A
        """
        self.encontrar_FEV()
        self.encontrar_entrante()
        self.encontrar_saliente()
        self.obtener_columna_pivote()
        self.modificar_linea_pivote()
        self.iteracion()

    def verificar_optimalidad(self):
        """ Verifica en la Fila U que no hayan negativos
            También revisa que no hayan soluciones múltiples
            E: N/A
            S: N/A
        """
        funcion_objetivo = self.matriz[1][1:]
        for valor in funcion_objetivo:
            if isinstance(valor, sympy.Basic) and valor != 0:
                funcion_objetivo[funcion_objetivo.index(valor)] = valor.subs({self.CONST_M:1000})
        #Primero revisa que no hayan negativos
        for n in funcion_objetivo[:-1]:
            if n < 0:
                return False

        #Si no los hay entonces revisa por soluciones múltiples      
        if (not(self.soluciones_multiples)):
            self.encontrar_basicas()
            i = 1
            while i < len(self.matriz[0])-2:
                if funcion_objetivo[i-1] == 0:
                    if not(self.matriz[0][i] in self.variables_basicas):
                        self.soluciones_multiples = True
                        return False
                i += 1
        
        return True

    def datos_solucion(self):
        """ Crea un string con los datos de solución de la matriz actual
            E: N/A
            S: string con los datos de la solución
        """
        str_matriz = "FEV: " + str(self.FEV)
        if self.is_max:
            str_matriz += "\nU: " + str(self.U)
        else:
            str_matriz += "\nU: " + str(self.U*-1)
        str_matriz += "\nVariable básica entrante: " + self.columna_pivote[0]
        str_matriz += "\nVariable básica saliente: " + self.fila_pivote[0]
        str_matriz += "\nPivote: " + str(self.pivote)

        return str_matriz

    def matriz_a_texto(self):
        """ Convierte la matriz en un string legible e imprimible
            E: N/A
            S: string de la matriz
        """

        matriz_redondeada = self.matriz
        fila = 1
        columna = 1
        while(fila < len(self.matriz)):
            columna = 1
            while(columna < len(self.matriz[0])):
                matriz_redondeada[fila][columna] = self.matriz[fila][columna]#.evalf(3)
                columna += 1
            fila += 1
        
        linea_string = [[str(casilla) for casilla in linea ] for linea in matriz_redondeada]
        lista_posicion = [max(map(len, columna)) for columna in zip(*linea_string)]
        cambia_formato = '\t\t'.join('{{:{}}}'.format(x) for x in lista_posicion)
        tabla = [cambia_formato.format(*linea) for linea in linea_string]

        return "\n".join(tabla)
        

    def datos_sol_optima(self):
        """ Crea un string con los datos de una solución óptima
            E: N/A
            S: string con datos de una solución óptima
        """
        self.encontrar_FEV()
        datos = "FEV: " + str(self.FEV)
        if self.is_max:
            datos += "\nU: " + str(self.U)
        else:
            datos += "\nU: " + str(self.U*-1)
        if self.degenerada:
            datos += "\nLa solución es degenerada"
        return datos

    def degeneradas_a_texto(self):
        i = 1
        datos = "Las variables "
        while i < len(self.var_degeneradas):
            if i == len(self.var_degeneradas)-1:
                datos += self.var_degeneradas[i]
            else:
                datos += self.var_degeneradas[i] + ", "
            i += 1
        datos += " son degeneradas, su resultado en común es " + str(round(self.var_degeneradas[0], 2))
        return datos

    def definir_ecuaciones_granm(self, diccionario_datos):
        
        self.CONST_M = Symbol('M')
        indice = 0
        diccionario_datos["fun_ob"] = [-x for x in diccionario_datos["fun_ob"]]
        var_basicas = []
        i = 0

        while i < len(diccionario_datos["simb_rest"]):
            if diccionario_datos["simb_rest"][i] == "<=":
                diccionario_datos["fun_ob"].append(Rational(0))
                var_basicas.append(len(diccionario_datos["fun_ob"]))

            elif diccionario_datos["simb_rest"][i] == ">=":
                diccionario_datos["fun_ob"].append(Rational(0))
                diccionario_datos["num_rest"] += 1 
                if diccionario_datos["optm"] == "max":
                    diccionario_datos["fun_ob"].append(self.CONST_M)
                    indice = diccionario_datos["fun_ob"].index(self.CONST_M, indice+1, len(diccionario_datos["fun_ob"]))
                    var_basicas.append(indice + 1)
                else:
                    diccionario_datos["fun_ob"].append(-self.CONST_M)
                    indice = diccionario_datos["fun_ob"].index(-self.CONST_M, indice+1, len(diccionario_datos["fun_ob"]))
                    var_basicas.append(indice + 1)
                self.var_artificiales.append(i)

            else:
                if diccionario_datos["optm"] == "max":
                    diccionario_datos["fun_ob"].append(self.CONST_M)
                    indice = diccionario_datos["fun_ob"].index(self.CONST_M, indice+1, len(diccionario_datos["fun_ob"]))
                    var_basicas.append(indice + 1)
                else:
                    diccionario_datos["fun_ob"].append(-self.CONST_M)
                    indice = diccionario_datos["fun_ob"].index(-self.CONST_M, indice+1, len(diccionario_datos["fun_ob"]))
                    var_basicas.append(indice + 1)
                self.var_artificiales.append(i)
            i += 1
        
        diccionario_datos["fun_ob"].append(Rational(0))

        if diccionario_datos["optm"] == "min":
            diccionario_datos["fun_ob"] = [-x for x in diccionario_datos["fun_ob"]]
            self.is_max = False

        
        for rest in diccionario_datos["rest"]: #Coloca los ceros en las restricciones
            tmp = rest.pop(-1)
            rest += [Rational(0)] * diccionario_datos["num_rest"]
            rest += [Rational(tmp)]
        
        i = diccionario_datos["num_var"]
        j = 0
        for rest in diccionario_datos["rest"]:
            if diccionario_datos["simb_rest"][j] == ">=":
                rest[i] = -1
                rest[i+1] = 1
                i += 1
            else:
                rest[i] = 1
            j += 1
            i += 1
        
        for i in self.var_artificiales:
            j = 0
            while j < len(diccionario_datos["fun_ob"]):
                diccionario_datos["fun_ob"][j] = diccionario_datos["fun_ob"][j] + (-self.CONST_M * diccionario_datos["rest"][i][j])
                j += 1

        self.crear_matriz([diccionario_datos["fun_ob"]] + diccionario_datos["rest"], var_basicas, len(diccionario_datos["fun_ob"])-1)

    def acomodar_diccionario(self, diccionario_datos):
        nueva_fun_ob = []
        nueva_simb_rest = []

        for fila in diccionario_datos["rest"]:
            nueva_fun_ob += [fila[-1]]

        nuevas_rest = []
        for fila in diccionario_datos["rest"]:
            nuevas_rest += [fila[:-1]]
        nuevas_rest = transpuesta(nuevas_rest)

        i = 0
        while i < len(nuevas_rest):
            nuevas_rest[i] += [diccionario_datos["fun_ob"][i]]
            i += 1

        if diccionario_datos["optm"] == "max":
            nueva_optm = "min"
            simbolo_restricciones = ">="
            nuevo_metodo = 1        #gran m
        else:
            nueva_optm = "max"
            simbolo_restricciones = "<="
            nuevo_metodo = 0        #simplex normal

        nuevo_num_var = diccionario_datos["num_rest"]

        nuevo_num_rest = diccionario_datos["num_var"]

        for i in nuevas_rest:
            nueva_simb_rest.append(simbolo_restricciones)
        
        #Guardamos los nuevos valores en el diccionario
        diccionario_datos["metodo"] = nuevo_metodo
        diccionario_datos["optm"] = nueva_optm
        diccionario_datos["num_var"] = nuevo_num_var
        diccionario_datos["num_rest"] = nuevo_num_rest
        diccionario_datos["fun_ob"] = nueva_fun_ob
        diccionario_datos["rest"] = nuevas_rest
        diccionario_datos["simb_rest"] = nueva_simb_rest
    
        return diccionario_datos


"""Fuera de clase Matriz"""

def leer_archivo(nombre_archivo):
    """ Función encargada de leer el archivo, también guarda todos los datos en un diccionario de datos
            E: string con el nombre del archivo
            S: el diccionario de datos con los datos ya guardados
    """
    contador = 0
    lista_datos = []
    diccionario_datos = { 
        "metodo" : 0, 
        "optm" : "", 
        "num_var" : 0, 
        "num_rest" : 0, 
        "fun_ob" : [], 
        "rest" : [],
        "simb_rest" : [] }
        
    try:
        with open(nombre_archivo,"r") as archivo:
            for lineas in archivo:
                lista_datos = lineas.split(",")
                
                if contador == 0:    
                    diccionario_datos["metodo"] = int(lista_datos[0])
                    diccionario_datos["optm"] = lista_datos[1]
                    diccionario_datos["num_var"] = int(lista_datos[2])
                    diccionario_datos["num_rest"] = int(lista_datos[3].replace("\n", ""))
                
                elif contador == 1:
                    diccionario_datos["fun_ob"] = list(map(Rational, lista_datos))  
     
                else:
                    diccionario_datos["simb_rest"].append(lista_datos[-2])
                    lista_datos.pop(-2)
                    diccionario_datos["rest"] += [list(map(Rational, lista_datos))]
                
                contador += 1
        archivo.close()
        return diccionario_datos
        
    except:
        print("\nEl archivo no se pudo abrir o no existe\n")
        quit()

def imprimir_ayuda():
    """ Imprime en pantalla la guia para la utilización del programa
            E: N/A
            S: N/A
    """
    str_ayuda = "\n   _____ _                 _           "
    str_ayuda += "\n  / ____(_)               | |          "
    str_ayuda += "\n | (___  _ _ __ ___  _ __ | | _____  __"
    str_ayuda += "\n  \___ \| | '_ ` _ \| '_ \| |/ _ \ \/ /"
    str_ayuda += "\n  ____) | | | | | | | |_) | |  __/>  < "
    str_ayuda += "\n |_____/|_|_| |_| |_| .__/|_|\___/_/\_\\"
    str_ayuda += "\n                    | |                "
    str_ayuda += "\n                    |_|\n"
    str_ayuda += "\nEjecución de programa: \n"
    str_ayuda += "\nPara ver la ayuda y correr el programa $ python simplex.py -h problema1.txt"
    str_ayuda += "\nPara ver la ayuda $ python simplex.py -h"
    str_ayuda += "\nPara solamente correr el programa $ python simplex.py problema1.txt\n"
    str_ayuda += "\nEstructura para formato en archivo de texto plano: \n"
    str_ayuda += "\nMétodo, optimización, número de variables de decisión, número de restricciones"
    str_ayuda += "\nCoeficientes de la función objetivo"
    str_ayuda += "\nCoeficientes de las restricciones y signo de restricción\n"
    str_ayuda += "\nMétodo es un valor numérico [ 0=Simplex, 1=GranM, 2=DosFases] "
    str_ayuda += "y optimización se indica de forma textual con min o max.\n"
    str_ayuda += "\nEjemplo para formato en archivo de texto plano: \n"
    str_ayuda += "\n---------------------------------------------------------------"
    str_ayuda += "\n| 0,max,2,3                                                   |"
    str_ayuda += "\n| 3,5                                                         |"
    str_ayuda += "\n| 2,1,<=,6                                                    |"
    str_ayuda += "\n| -1,3,<=,9                                                   |"
    str_ayuda += "\n| 0,1,<=,4                                                    |"
    str_ayuda += "\n---------------------------------------------------------------\n"

    print(str_ayuda)

def escribir_archivo(nombre_archivo, texto):
    """ Función encargada de escribir texto en un archivo
            E: recibe la ruta del archivo y el texto a escribir
            S: N/A
    """
    nombre_archivo = str(nombre_archivo).replace(".txt", "")
    nombre_archivo += "_solution.txt"
    try:
        with open(nombre_archivo,"a") as archivo:
            archivo.write(texto + os.linesep)

    except:
        print("\nNo se pudo crear o abrir el archivo\n")
    
    archivo.close()

def limpiar_archivo_solucion(nombre_archivo):
    """ Limpia el archivo en caso de que tenga algo escrito
            E: ruta del archivo
            S: N/A
    """
    nombre_archivo= str(nombre_archivo).replace(".txt", "")
    nombre_archivo += "_solution.txt"
    try:
        with open(nombre_archivo,"w") as archivo:
            archivo.write("")

        archivo.close()
    except:
        print("\nNo se pudo crear o abrir el archivo\n")

def obtener_solucion(nombre_archivo):
    num_iteracion = 0
    diccionario_datos = leer_archivo(nombre_archivo)
    matriz = Matriz(diccionario_datos)
    limpiar_archivo_solucion(nombre_archivo)
    matriz.nom_archivo = nombre_archivo
    
    while(True):
        
        if matriz.soluciones_multiples:
            print(matriz.datos_sol_optima())

        escribir_archivo(nombre_archivo,"\nIteracion " + str(num_iteracion))
        escribir_archivo(nombre_archivo,matriz.matriz_a_texto())
        matriz.iterar()
        escribir_archivo(nombre_archivo,matriz.datos_solucion())

        if (matriz.verificar_optimalidad()):
            if matriz.soluciones_multiples:
                escribir_archivo(nombre_archivo, "Solución múltiple en: " + matriz.columna_pivote[0])
                print ("Solución múltiple en: " + matriz.columna_pivote[0])
                escribir_archivo(nombre_archivo,"\nIteracion extra")
                print("\nIteracion extra")
            else:
                escribir_archivo(nombre_archivo,"\nIteracion Final")  
            escribir_archivo(nombre_archivo, matriz.matriz_a_texto())
            print(matriz.datos_sol_optima())
            escribir_archivo(nombre_archivo,matriz.datos_sol_optima())
            break

        num_iteracion += 1

def transpuesta(matriz):
    matriz_transpuesta = []

    for j in range(len(matriz[0])):
        matriz_transpuesta.append([])
        for i in range(len(matriz)):
            matriz_transpuesta[j].append(matriz[i][j])

    return matriz_transpuesta

def principal(args):
    """ Función encargada de la ejecución del programa
            E: recibe los argumento ingresados por consola
            S: N/A
    """

    if len(args) == 3 and args[1] == "-h":
        imprimir_ayuda()
        obtener_solucion(args[2])

    elif len(args) == 2:

        if args[1] == "-h":
            imprimir_ayuda()
        
        else:
            obtener_solucion(args[1])
    else:
        print("\nIngrese [-h] para recibir ayuda de utilización del programa\n")
        return

principal(sys.argv)