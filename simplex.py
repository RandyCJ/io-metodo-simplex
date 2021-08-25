import sys
from iteration_utilities import duplicates

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

def crear_matriz(matriz, variables_basicas, cant_variables):
    
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
    nueva_matriz = [encabezado] + nueva_matriz
    return nueva_matriz

def definir_ecuaciones(diccionario_datos):
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
    return crear_matriz([diccionario_datos["fun_ob"]] + diccionario_datos["rest"], var_basicas, diccionario_datos["num_var"] + diccionario_datos["num_rest"])

def encontrar_entrante(matriz):
    fila = matriz[1][1:]
    fila.sort()
    columna = matriz[1].index(fila[0])
    entrante = matriz[0][matriz[1].index(fila[0])]
    return (entrante, columna)


def encontrar_saliente(matriz, entrante):
    saliente = ("",0)
    lista = []
    lista_degenerados = []
    acotada = True

    for i in matriz[2:]:
        if i[entrante[1]] > 0 and i[-1] > 0:
            acotada = False
            lista += [(i[-1]/i[entrante[1]],i[0])]
            lista_degenerados += [i[-1]/i[entrante[1]]]

    if acotada == True:
        saliente = "Es acotada"

    else:
        lista.sort()
        saliente = lista[0][1]

        ind_saliente = -1
        for fila in matriz:
            if fila[0] != saliente:
                ind_saliente += 1
        saliente = (lista[0][1],ind_saliente)

    if list(duplicates(lista_degenerados)) != []:
        saliente = "Es degenerada"

    return saliente

def encontrar_FEV(matriz, diccionario_datos):
    U = matriz[1][-1]
    columna = 1
    fila = 1
    lista_FEV = []

    for i in matriz[1][1:]:
        if i == 0:
            fila = 1
            while (columna < diccionario_datos["num_rest"]+4 and fila < diccionario_datos["num_rest"]+2):
                if matriz[fila][columna] == 1:
                    lista_FEV += [matriz[fila][-1]]   

                fila += 1
        else:
            lista_FEV += [0]
        columna += 1

    return (U, lista_FEV)

def llenar_columna(matriz, entrante):
    i = 1
    while i < len(matriz):
        matriz[i][entrante[1]] = 0
        i += 1

    return matriz

def llenar_fila(pivote, entrante, saliente, nueva_matriz):
    i = 1
    nueva_matriz[saliente[1]][0] = entrante[0]
    while i < len(nueva_matriz[0]):
        nueva_matriz[saliente[1]][i] = nueva_matriz[saliente[1]][i]/pivote
        i += 1
    nueva_matriz[saliente[1]][entrante[1]] = 1
    return nueva_matriz

def principal(args):
    
    diccionario_datos = {}

    if len(args) == 3 and args[1] == "-h":
        print ("\nExiste el argumento de ayuda y el argumento del archivo\n")
        diccionario_datos = leer_archivo(args[2])
        matriz = definir_ecuaciones(diccionario_datos)
        FEV = encontrar_FEV(matriz, diccionario_datos)
        entrante = encontrar_entrante(matriz)
        saliente = encontrar_saliente(matriz, entrante)
        pivote = matriz[saliente[1]][entrante[1]]
        nueva_matriz = llenar_columna(matriz, entrante)
        nueva_matriz = llenar_fila(pivote, entrante, saliente, nueva_matriz)
        for fila in nueva_matriz:
            print(fila)
        
    
    elif len(args) == 2:

        if args[1] == "-h":
            print("\nAca escribiremos el código de ayuda\n")
        
        else:
            print("\nSolo el argumento del archivo\n")
            diccionario_datos = leer_archivo(args[1])
            matriz = definir_ecuaciones(diccionario_datos)
            entrante = encontrar_entrante(matriz)

    else:
        print("\nIngrese [-h] para recibir ayuda de utilización del programa\n")
        return


principal(sys.argv)