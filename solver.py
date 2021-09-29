import sys 
import os
import sympy
from iteration_utilities import duplicates
from sympy import *
from sympy.core.rules import Transform
from simplex import Matriz

def crear_matriz(matriz, variables_basicas, cant_variables, es_maximizacion):
    """ Crea la matriz junto con el encabezado de las variables y los números y valores de cada una
        E: una matriz con solo los valores, las variables que van en la columna 0 y la cantidad de variables que deben haber, un booleano si es Max o Min
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
            if es_maximizacion:
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

    matriz = [encabezado] + nueva_matriz
    return matriz

def definir_ecuaciones(diccionario_datos, CONST_M):
    """ Convierte el diccionario de datos a las ecuaciones que se utilizaran para las tablas
        E: Recibe el diccionario de datos con los datos que se recolectaron al leer el archivo
        S: N/A
    """
    dual = False

    if diccionario_datos["metodo"] == 3:   #dual
        dual = True
        diccionario_datos = acomodar_diccionario(diccionario_datos)

    diccionario_datos["fun_ob"] = [-x for x in diccionario_datos["fun_ob"]] #Se vuelven negativos todos los números de la función objetivo

    if diccionario_datos["metodo"] == 1:   #gran m
        return (definir_ecuaciones_granm(diccionario_datos, CONST_M),dual)
    
    if diccionario_datos["metodo"] == 2:   #Dos fases
        return (definir_ecuaciones_primera_fase(diccionario_datos),dual)

    diccionario_datos["fun_ob"] += [Rational(0)] * (diccionario_datos["num_rest"] + 1) #Agrega los ceros dependiendo de la cantidad de restricciones

    for rest in diccionario_datos["rest"]: #Coloca los ceros y unos en las restricciones
        tmp = rest.pop(-1)
        rest += [Rational(0)] * diccionario_datos["num_rest"]
        rest += [Rational(tmp)]

    i = diccionario_datos["num_var"]
    for rest in diccionario_datos["rest"]:
        rest[i] = Rational(1)
        i += 1

    var_basicas = [x for x in range(diccionario_datos["num_var"]+1, len(diccionario_datos["fun_ob"]))]
    return ([crear_matriz([diccionario_datos["fun_ob"]] + diccionario_datos["rest"], var_basicas, diccionario_datos["num_var"] + diccionario_datos["num_rest"], True)],dual)

def definir_ecuaciones_granm(diccionario_datos, CONST_M):
    indice = 0
    var_basicas = []
    indice_var_artificiales = []
    var_artificiales = []
    i = 0
    es_maximizacion = True

    while i < len(diccionario_datos["simb_rest"]):
        if diccionario_datos["simb_rest"][i] == "<=":
            diccionario_datos["fun_ob"].append(Rational(0))
            var_basicas.append(len(diccionario_datos["fun_ob"]))

        elif diccionario_datos["simb_rest"][i] == ">=":
            diccionario_datos["fun_ob"].append(Rational(0))
            diccionario_datos["num_rest"] += 1 
            if diccionario_datos["optm"] == "max":
                diccionario_datos["fun_ob"].append(CONST_M)
                indice = diccionario_datos["fun_ob"].index(CONST_M, indice+1, len(diccionario_datos["fun_ob"]))
            else:
                diccionario_datos["fun_ob"].append(-CONST_M)
                indice = diccionario_datos["fun_ob"].index(-CONST_M, indice+1, len(diccionario_datos["fun_ob"]))
            var_basicas.append(indice + 1)
            indice_var_artificiales.append(i)
            var_artificiales.append(indice + 1)

        else:
            if diccionario_datos["optm"] == "max":
                diccionario_datos["fun_ob"].append(CONST_M)
                indice = diccionario_datos["fun_ob"].index(CONST_M, indice+1, len(diccionario_datos["fun_ob"]))
            else:
                diccionario_datos["fun_ob"].append(-CONST_M)
                indice = diccionario_datos["fun_ob"].index(-CONST_M, indice+1, len(diccionario_datos["fun_ob"]))
            var_basicas.append(indice + 1)
            indice_var_artificiales.append(i)
            var_artificiales.append(indice + 1)
        i += 1
    
    diccionario_datos["fun_ob"].append(Rational(0))

    if diccionario_datos["optm"] == "min":
        diccionario_datos["fun_ob"] = [-x for x in diccionario_datos["fun_ob"]]
        es_maximizacion = False
    
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
    
    for i in indice_var_artificiales:
        j = 0
        while j < len(diccionario_datos["fun_ob"]):
            diccionario_datos["fun_ob"][j] = diccionario_datos["fun_ob"][j] + (-CONST_M * diccionario_datos["rest"][i][j])
            j += 1
    return [crear_matriz([diccionario_datos["fun_ob"]] + diccionario_datos["rest"], var_basicas, len(diccionario_datos["fun_ob"])-1, es_maximizacion), var_artificiales]
def definir_ecuaciones_primera_fase(diccionario_datos):     
    var_basicas = []
    var_artificiales = []
    var_exceso = [] 
    fun_ob_1 = []
    rest_artificiales = []
    resultado_rest = 0
    i = 0
    indice = 0
    agregar_0 = 0
    variable = 0
    es_maximizacion = True
    fun_ob_1.extend([0] * len (diccionario_datos["rest"][0]))
    
    if diccionario_datos["optm"] == "min":
        es_maximizacion = False

    while i < len(diccionario_datos["simb_rest"]):

        if diccionario_datos["simb_rest"][i] == "<=":
            
            cantidad = len(var_exceso) + len (var_artificiales)

            while agregar_0 < cantidad :
                resultado_rest = diccionario_datos["rest"][i][-1]
                diccionario_datos["rest"][i].pop(-1)
                diccionario_datos["rest"][i].append(0)
                diccionario_datos["rest"][i].append(resultado_rest)
                agregar_0 += 1 
            
            resultado_rest = diccionario_datos["rest"][i][-1]
            diccionario_datos["rest"][i].pop(-1)
            diccionario_datos["rest"][i].append(1)
            var_exceso.append(len(diccionario_datos["rest"][i])-1)
            var_basicas.append(len(diccionario_datos["rest"][i])-1)
            
            diccionario_datos["rest"][i].append(resultado_rest)
            
            if  i-1 >= 0:
                for restriccion in diccionario_datos["rest"]:
                    if len (diccionario_datos["rest"][i]) > len (diccionario_datos["rest"][i-1]):
                        extra_ceros = len (diccionario_datos["rest"][i]) - len (diccionario_datos["rest"][i-1])
                        resultado_rest= restriccion.pop(-1)
                        restriccion.extend([0] * extra_ceros)
                        restriccion.append(resultado_rest)
            
            fun_ob_1.extend([0])
            agregar_0 = 0
            
        
        elif diccionario_datos["simb_rest"][i] == ">=":
            
            cantidad = len(var_exceso) + len (var_artificiales)

            while agregar_0 < cantidad :
                resultado_rest = diccionario_datos["rest"][i][-1]
                diccionario_datos["rest"][i].pop(-1)
                diccionario_datos["rest"][i].append(0)
                diccionario_datos["rest"][i].append(resultado_rest)
                agregar_0 += 1  
            
            
            resultado_rest = diccionario_datos["rest"][i][-1]
            diccionario_datos["rest"][i].pop(-1)
            diccionario_datos["rest"][i].append(-1)
            var_exceso.append(len(diccionario_datos["rest"][i])-1)
            diccionario_datos["rest"][i].append(1)
            var_artificiales.append(len(diccionario_datos["rest"][i])-1)
            var_basicas.append(len(diccionario_datos["rest"][i])-1)

            diccionario_datos["rest"][i].append(resultado_rest)

            if  i-1 >= 0:
                for restriccion in diccionario_datos["rest"]:
                    if len (diccionario_datos["rest"][i]) > len (diccionario_datos["rest"][i-1]):
                        extra_ceros = len (diccionario_datos["rest"][i]) - len (diccionario_datos["rest"][i-1])
                        resultado_rest= restriccion.pop(-1)
                        restriccion.extend([0] * extra_ceros)
                        restriccion.append(resultado_rest)
            
            
            resultado_rest = fun_ob_1.pop(-1) + diccionario_datos["rest"][i][-1]
            fun_ob_1.extend([0] * 2)
        
            while variable <  len(diccionario_datos["rest"][i])-1:
                if variable in var_artificiales:
                    fun_ob_1[variable] = 0
                
                elif variable in var_exceso:
                    #print(variable)
                    fun_ob_1[variable] = 1

                elif diccionario_datos["rest"][i][variable] >= 0:
                    fun_ob_1[variable] = fun_ob_1[variable] - (diccionario_datos["rest"][i][variable])    
                
                elif diccionario_datos["rest"][i][variable] < 0:
                    fun_ob_1[variable] = fun_ob_1[variable] + (diccionario_datos["rest"][i][variable]) 
                
                variable += 1
            fun_ob_1.append(resultado_rest)
            agregar_0 = 0
            variable = 0
            rest_artificiales.append(i)
            
            
        elif diccionario_datos["simb_rest"][i] == "=":
            
            cantidad = len(var_exceso) + len (var_artificiales)

            while agregar_0 < cantidad :
                resultado_rest = diccionario_datos["rest"][i][-1]
                diccionario_datos["rest"][i].pop(-1)
                diccionario_datos["rest"][i].append(0)
                diccionario_datos["rest"][i].append(resultado_rest)
                agregar_0 += 1  
            
            
            resultado_rest = diccionario_datos["rest"][i][-1]
            diccionario_datos["rest"][i].pop(-1)
            diccionario_datos["rest"][i].append(1)
            var_artificiales.append(len(diccionario_datos["rest"][i])-1)
            var_basicas.append(len(diccionario_datos["rest"][i])-1)
            diccionario_datos["rest"][i].append(resultado_rest)

            if  i-1 >= 0:
                for restriccion in diccionario_datos["rest"]:
                    if len (diccionario_datos["rest"][i]) > len (diccionario_datos["rest"][i-1]):
                        extra_ceros = len (diccionario_datos["rest"][i]) - len (diccionario_datos["rest"][i-1])
                        resultado_rest= restriccion.pop(-1)
                        restriccion.extend([0] * extra_ceros)
                        restriccion.append(resultado_rest)
                
                

            resultado_rest = fun_ob_1.pop(-1) + diccionario_datos["rest"][i][-1]
            fun_ob_1.extend([0])
        
            while variable <  len(diccionario_datos["rest"][i])-1:
                if variable in var_artificiales:
                    fun_ob_1[variable] = 0
                
                elif variable in var_exceso:
                    fun_ob_1[variable] = 1

                elif diccionario_datos["rest"][i][variable] >= 0:
                    fun_ob_1[variable] = fun_ob_1[variable] - (diccionario_datos["rest"][i][variable])    
                
                elif diccionario_datos["rest"][i][variable] < 0:
                    fun_ob_1[variable] = fun_ob_1[variable] + (diccionario_datos["rest"][i][variable]) 
            
                variable += 1
            fun_ob_1.append(resultado_rest)
            agregar_0 = 0
            variable = 0
            rest_artificiales.append(i)
        
        i+=1
    fun_ob_1[-1]= -fun_ob_1[-1]

    while indice < len(var_basicas):
        var_basicas[indice] += 1 
        indice+=1
    indice = 0
    while indice < len(var_artificiales):
        var_artificiales[indice] += 1 
        indice+=1

    return [crear_matriz([fun_ob_1] + diccionario_datos["rest"], var_basicas , len(fun_ob_1)-1, es_maximizacion), var_artificiales]


def acomodar_diccionario(diccionario_datos):
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

def manejar_no_acotada(matriz, nombre_archivo):
    if matriz.dual:
        msj_acotada = 'La solución dual es no acotada en '
        msj_acotada += 'la columna ' + str(matriz.columna_pivote[1]+1)
        msj_acotada += ' con la VB entrante ' + matriz.columna_pivote[0]
        msj_acotada += '\nPor ello la solución primal no tiene soluciones factibles'
    else:
        msj_acotada = 'La columna ' + str(matriz.columna_pivote[1]+1)
        msj_acotada += ' con la VB entrante ' + matriz.columna_pivote[0]
        msj_acotada += ' no tiene números mayores a 0'
        msj_acotada += "\nPor ello la solución es no acotada"

    escribir_archivo(nombre_archivo,"\n" + msj_acotada)
    print(msj_acotada)
    quit()

def manejar_no_factible(matriz, nombre_archivo, var_no_factible):
    if matriz.dual:
        msj_no_factible = "La variable artificial " + var_no_factible + " es positiva "
        msj_no_factible += "esto hace la solución dual no factible"
        msj_no_factible += "\nPor ello la solución primal no tiene soluciones factibles "
        msj_no_factible += "o es no acotada"
    else:
        msj_no_factible = "La variable artificial " + var_no_factible + " es positiva\n"
        msj_no_factible += "Por lo tanto la solución no es factible"
    escribir_archivo(nombre_archivo, "\n" + msj_no_factible)
    print(msj_no_factible)

def realizar_iteraciones(matriz, nombre_archivo):

    num_iteracion = 0

    while(True):
        
        if matriz.soluciones_multiples and not(matriz.dual):
            print(matriz.datos_sol_optima())
        
        escribir_archivo(nombre_archivo,"\nIteracion " + str(num_iteracion))
        escribir_archivo(nombre_archivo,matriz.matriz_a_texto())
        try:
            matriz.iterar()
        except Exception as e:
            manejar_no_acotada(matriz, nombre_archivo)
        escribir_archivo(nombre_archivo,matriz.datos_solucion())
            
        if matriz.verificar_optimalidad():

            if matriz.dual:
                escribir_archivo(nombre_archivo,"\nIteracion Final Dual")
                escribir_archivo(nombre_archivo, matriz.matriz_a_texto())
                print("Solución dual:")
                print(matriz.datos_sol_optima())
                print("\nSolución primal:")
                print(datos_sol_optima_dual(matriz))
                escribir_archivo(nombre_archivo,"\nSolución dual:")
                escribir_archivo(nombre_archivo,matriz.datos_sol_optima())
                escribir_archivo(nombre_archivo,"\nSolución primal:")
                escribir_archivo(nombre_archivo,datos_sol_optima_dual(matriz))
                no_factible = matriz.verificar_artificiales()
                if no_factible != 0:
                   manejar_no_factible(matriz,nombre_archivo, no_factible)
                break

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
            no_factible = matriz.verificar_artificiales()
            if no_factible != 0:
                manejar_no_factible(matriz, nombre_archivo, no_factible)
            break

        num_iteracion += 1

    return matriz

def obtener_solucion(nombre_archivo):
    CONST_M = Symbol('M')
    diccionario_datos = leer_archivo(nombre_archivo)
    tupla_tmp = definir_ecuaciones(diccionario_datos, CONST_M)
    matriz_inicial = tupla_tmp[0]
    matriz = Matriz(matriz_inicial[0])
    matriz.dual = tupla_tmp[1]

    if diccionario_datos["metodo"] == 1:
        matriz.definir_artificiales(matriz_inicial[1], CONST_M)

    if diccionario_datos["metodo"] == 2:
        matriz.definir_artificiales(matriz_inicial[1], CONST_M)
        matriz.dos_fases = True
    
    sacar_holgura(matriz, diccionario_datos)
    limpiar_archivo_solucion(nombre_archivo)
    matriz.nom_archivo = nombre_archivo
    
    matriz = realizar_iteraciones(matriz, nombre_archivo)
    
    if matriz.dos_fases:
        matriz.matriz = eliminar_artificiales(matriz.matriz,matriz.var_artificiales)
        matriz.matriz = cambiar_fun_obj(matriz,diccionario_datos,matriz.es_max)
        matriz = realizar_iteraciones(matriz, nombre_archivo)

def eliminar_artificiales(matriz, var_artificiales):
    i=1
    nueva_matriz = transpuesta(matriz)
    
    while i < len(nueva_matriz)-1:
        if nueva_matriz[i][0] in var_artificiales:
            nueva_matriz.pop(i)
            i-=1
        i+=1

    nueva_matriz = transpuesta(nueva_matriz)
    return nueva_matriz

def cambiar_fun_obj(obj_matriz, diccionario_datos, es_maximizacion):
    matriz = obj_matriz.matriz
    matriz = transpuesta(matriz)
    var_basicas = matriz[0][2:]
    matriz = transpuesta(matriz)
    fun_ob_2 = []
    
    if es_maximizacion:
        fun_ob_2 = ["U"] + diccionario_datos["fun_ob"]
        matriz[1] = fun_ob_2
    else:
        i=0
        while i < len(diccionario_datos["fun_ob"]):
            fun_ob_2.append(-diccionario_datos["fun_ob"][i])
            i+=1
        fun_ob_2 = ["-U"] + fun_ob_2
    
    cantidad_0 = len(matriz[0]) - len (fun_ob_2)
    fun_ob_2.extend([0]*cantidad_0)
    matriz[1]= fun_ob_2

    if obj_matriz.verificar_optimalidad():
        for variable in matriz[0]:
            
            if variable in var_basicas:
                h = 0
                while h < len(matriz):    
                    if variable == matriz[h][0]:
                        fila_multiplicacion = matriz[h]
                
                        realizar_0 = matriz[1][matriz[0].index(variable)]
                        i = 1

                        for numero in matriz[1][1:]:

                            matriz[1][i] = numero - (realizar_0 * fila_multiplicacion[i])
                            i+=1
                    h+=1
    return matriz


def transpuesta(matriz):
    matriz_transpuesta = []

    for j in range(len(matriz[0])):
        matriz_transpuesta.append([])
        for i in range(len(matriz)):
            matriz_transpuesta[j].append(matriz[i][j])

    return matriz_transpuesta

def encontrar_FEV_dual(matriz):

        matriz.U = matriz.matriz[1][-1]
        matriz.FEV = []
        columna = 1
        tmp = []                        

        while columna < len(matriz.matriz[0][:-1]):
            if matriz.matriz[0][columna] in matriz.var_holgura:
                tmp.append(matriz.matriz[1][columna])
            columna += 1
        matriz.FEV += tmp

def datos_sol_optima_dual(matriz):

        encontrar_FEV_dual(matriz)
        i = 0
        datos = "BF: \n"
        while i < len(matriz.FEV):
            datos += "X" + str(i+1) + ": " + str(matriz.FEV[i]) + " "
            i += 1
        if matriz.es_max:
            datos += "\nU: " + str(matriz.U)
        else:
            datos += "\nU: " + str(matriz.U*-1)

        return datos

def sacar_holgura(matriz, diccionario_datos):

    i = diccionario_datos["num_var"]
    candidatos_holgura = []
    while i < len(matriz.matriz[0][1:-1]):
        candidatos_holgura.append("X" + str(i+1))
        i += 1

    i = 0
    while i < len(matriz.matriz[0][1:]):
        if matriz.matriz[0][i] in candidatos_holgura:
            matriz.var_holgura.append(matriz.matriz[0][i])
        i += 1

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