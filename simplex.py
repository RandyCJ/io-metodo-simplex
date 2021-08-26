import sys
from iteration_utilities import duplicates

class Matriz:
    
    matriz = []
    variables_basicas = []
    columna_pivote = ()
    fila_pivote = ()
    indice_pivote = ()
    pivote = 1
    acotada = False
    degenerada = False
    var_degeneradas = []
    soluciones_multiples = False
    U = 0
    FEV = []

    def __init__(self, dict_datos) -> None:
        self.definir_ecuaciones(dict_datos)

    def crear_matriz(self, matriz, variables_basicas, cant_variables):
        encabezado = ["VB"]

        for i in range(cant_variables):
            encabezado.append("X" + str(i+1))
        encabezado.append("LD")
        
        i = 0
        nueva_matriz = []

        for fila in matriz:
            nueva_fila = []
            if i == 0:
                nueva_fila.append("U")
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
        """ Convierte el diccionario de datos a las ecuaciones que se utilizaran para las tablas"""

        diccionario_datos["fun_ob"] = [-x for x in diccionario_datos["fun_ob"]] #Se vuelven negativos todos los números
        diccionario_datos["fun_ob"] += [0] * (diccionario_datos["num_rest"] + 1)

        for rest in diccionario_datos["rest"]:
            tmp = rest.pop(-1)
            rest += [0] * diccionario_datos["num_rest"]
            rest += [tmp]

        i = diccionario_datos["num_var"]
        for rest in diccionario_datos["rest"]:
            rest[i] = 1
            i += 1

        var_basicas = [x for x in range(diccionario_datos["num_var"]+1, (len(diccionario_datos["fun_ob"])-1)+1)]
        self.crear_matriz([diccionario_datos["fun_ob"]] + diccionario_datos["rest"], var_basicas, diccionario_datos["num_var"] + diccionario_datos["num_rest"])

    def encontrar_entrante(self):
        if self.soluciones_multiples: #Si es multiple no hace falta sacar columna entrante
            return 0
        fila = self.matriz[1][1:]
        fila.sort()
        self.columna_pivote = (self.matriz[0][self.matriz[1].index(fila[0])], self.matriz[1].index(fila[0]))

    def encontrar_saliente(self):
        lista_divisiones = []
        lista_degenerados = []
        self.var_degeneradas = []
        self.acotada = True
        vb_saliente = ""

        for i in self.matriz[2:]:
            if i[self.columna_pivote[1]] > 0 and i[-1] >= 0:
                self.acotada = False
                lista_divisiones += [(i[-1]/i[self.columna_pivote[1]], i[0])]
                lista_degenerados += [i[-1]/i[self.columna_pivote[1]]]

        if self.acotada:
            print (f'La columna {self.columna_pivote[1]+1!r} ' +
                   f'con la VB entrante {self.columna_pivote[0]} no tiene números mayores a 0')
            print("Por ello la solución es no acotada")
            quit()

        lista_divisiones.sort()
        vb_saliente = lista_divisiones[0][1]
        ind_saliente = 0

        for fila in self.matriz:
            if fila[0] != vb_saliente:
                ind_saliente += 1
                continue
            break

        self.fila_pivote = (vb_saliente, ind_saliente)
         
        if lista_divisiones[0][0] in list(duplicates(lista_degenerados)):
            self.degenerada = True
            i = 0
            self.var_degeneradas = [lista_divisiones[0][0]]
            while i < len(lista_divisiones):
                if lista_divisiones[i][0] == lista_divisiones[0][0]:
                    self.var_degeneradas += [lista_divisiones[i][1]]
                i += 1
    
    def encontrar_basicas(self):
        var_basicas = []
        for fila in self.matriz[2:]:
            var_basicas.append(fila[0])
        self.variables_basicas = var_basicas

    def encontrar_FEV(self):
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
        """Almacena la columna pivote donde está el pivote"""
        fila = 1
        columna = []
        while (fila < len(self.matriz)):
            columna += [-self.matriz[fila][self.columna_pivote[1]]]
            fila += 1
        self.columna_pivote = (self.columna_pivote[0], self.columna_pivote[1], columna)

    def modificar_linea_pivote(self):
        #Primero se convierten todos los números de la columna donde está el pivote a 0
        self.pivote = self.matriz[self.fila_pivote[1]][self.columna_pivote[1]]
        i = 1
        while i < len(self.matriz):
            self.matriz[i][self.columna_pivote[1]] = 0
            i += 1

        #Después se divide toda la linea entre el pivote
        i = 1
        self.matriz[self.fila_pivote[1]][0] = self.columna_pivote[0] #Sale VB Saliente y entra la VB Entrante
        while i < len(self.matriz[0]):
            self.matriz[self.fila_pivote[1]][i] = self.matriz[self.fila_pivote[1]][i]/self.pivote
            i += 1
        self.matriz[self.fila_pivote[1]][self.columna_pivote[1]] = 1

    def iteracion(self):
        num_fila = 0
        for fila in self.matriz:
            if num_fila != self.fila_pivote[1] and num_fila != 0 :
                indice = 0
                for num in fila[1:]:
                    if indice != self.columna_pivote[1] - 1 and indice <= len(self.matriz) + 1:
                        self.matriz[num_fila][indice + 1] = float(num + self.columna_pivote[2][num_fila - 1] * self.matriz[self.fila_pivote[1]][indice + 1])
                    indice += 1                
            num_fila += 1
    
    def iterar(self):
        self.encontrar_FEV()
        self.encontrar_entrante()
        self.encontrar_saliente()
        self.obtener_columna_pivote()
        self.modificar_linea_pivote()
        self.iteracion()

    def verificar_optimalidad(self):
        funcion_objetivo = self.matriz[1]
        #Primero revisa que no hayan negativos
        for n in funcion_objetivo[1:]:
            if n < 0:
                return False

        #Si no los hay entonces revisa por soluciones múltiples      
        if (not(self.soluciones_multiples)):
            self.encontrar_basicas()
            i = 1
            while i < len(self.matriz[0])-1:
                if funcion_objetivo[i] == 0:
                    if not(self.matriz[0][i] in self.variables_basicas):
                        self.soluciones_multiples = True
                        print("SM en " + str(self.matriz[0][i]))
                        return False
                i += 1
            
        return True

    def toString(self):

        str_matriz = "\nFEV: " + str(self.FEV)
        str_matriz += "\nU: " + str(self.U)
        str_matriz += "\nVariable básica entrante: " + self.columna_pivote[0]
        str_matriz += "\nVariable básica saliente: " + self.fila_pivote[0]
        str_matriz += "\nPivote: " + str(self.pivote)

        return str_matriz

    def imprimir_matriz(self):
        print("\n")
        for i in self.matriz:
            print(i)
    
    def datos_sol_optima(self):
        self.encontrar_FEV()
        datos = "\nFEV: " + str(self.FEV)
        datos += "\nU: " + str(self.U)
        if self.degenerada:
            datos += "\nLa solución es degenerada"
        return datos

"""Fuera de clase Matriz"""

def leer_archivo(nombre_archivo):
    
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
                    diccionario_datos["fun_ob"] = list(map(int, lista_datos))  
     
                else:
                    diccionario_datos["simb_rest"].append(lista_datos[-2])
                    lista_datos.pop(-2)
                    diccionario_datos["rest"] += [list(map(int, lista_datos))]
                
                contador += 1

        return diccionario_datos
        
    except:
        print("\nEl archivo no se pudo abrir o no existe\n")

def imprimir_ayuda():
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
    str_ayuda += "\nEjemplo para formato de archivo de datos: \n"
    str_ayuda += "\n---------------------------------------------------------------"
    str_ayuda += "\n| 0,max,2,3                                                   |"
    str_ayuda += "\n| 3,5                                                         |"
    str_ayuda += "\n| 2,1,<=,6                                                    |"
    str_ayuda += "\n| -1,3,<=,9                                                   |"
    str_ayuda += "\n| 0,1,<=,4                                                    |"
    str_ayuda += "\n---------------------------------------------------------------\n"
    str_ayuda += "\nSe debe utilizar un archivo de texto plano\n"

    print(str_ayuda)

def principal(args):
    
    diccionario_datos = {}

    if len(args) == 3 and args[1] == "-h":
        imprimir_ayuda()
        diccionario_datos = leer_archivo(args[2])
        matriz = Matriz(diccionario_datos)

        while(True):
            matriz.imprimir_matriz()
            matriz.iterar()
            print(matriz.toString())
            
            if (matriz.verificar_optimalidad()):
                if matriz.soluciones_multiples:
                    print("\nIteracion extra")
                matriz.imprimir_matriz()
                print(matriz.datos_sol_optima())
                break

    elif len(args) == 2:

        if args[1] == "-h":
            imprimir_ayuda()
        
        else:
            diccionario_datos = leer_archivo(args[1])
            matriz = Matriz(diccionario_datos)

            while(True):
                matriz.imprimir_matriz()
                matriz.iterar()
                print(matriz.toString())
                
                if (matriz.verificar_optimalidad()):
                    if matriz.soluciones_multiples:
                        print("\nIteracion extra")
                    matriz.imprimir_matriz()
                    print(matriz.datos_sol_optima())
                    break
    else:
        print("\nIngrese [-h] para recibir ayuda de utilización del programa\n")
        return

principal(sys.argv)