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
    lista_divisiones = []
    lista_degenerados = []
    acotada = True

    for i in matriz[2:]:
        if i[entrante[1]] > 0 and i[-1] > 0: #No sabemos si se puede dividir 0/n
            acotada = False
            lista_divisiones += [(i[-1]/i[entrante[1]], i[0])]
            lista_degenerados += [i[-1]/i[entrante[1]]]

    if acotada:
        saliente = "Es acotada"

    else:
        lista_divisiones.sort()
        saliente = lista_divisiones[0][1]
        ind_saliente = 0
        for fila in matriz:
            if fila[0] != saliente:
                ind_saliente += 1
                continue
            break

        saliente = (saliente, ind_saliente)

    
    if lista_divisiones[0][0] in list(duplicates(lista_degenerados)):
        saliente = "Es degenerada"
        print("La solución es degenerada")
        quit()

    return saliente

def encontrar_FEV(matriz):
    U = matriz[1][-1]
    columna = 1
    lista_FEV = []
    var_basicas = encontrar_basicas(matriz)
    lista_multiples = []
    bandera_multiples = False

    for i in matriz[1][1:]:
        if i == 0:
            if not(matriz[0][columna] in var_basicas) and columna < len(matriz[0])-1:
                bandera_multiples = True
                lista_multiples = [bandera_multiples, matriz[0][columna], columna]
                lista_FEV += [0]
            else:
                fila = 1
                while (columna < len(matriz[0]) and fila < len(matriz)):
                    if matriz[fila][columna] == 1:
                        lista_FEV += [matriz[fila][-1]]   
                    
                    fila += 1

        elif columna < len(matriz[0])-1:
            lista_FEV += [0]

        columna += 1
        
    return (U, lista_FEV, lista_multiples)

def llenar_columna(matriz, entrante):
    nueva_matriz = matriz
    i = 1
    while i < len(nueva_matriz):
        nueva_matriz[i][entrante[1]] = 0
        i += 1

    return nueva_matriz

def llenar_fila(pivote, entrante, saliente, nueva_matriz):
    i = 1
    nueva_matriz[saliente[1]][0] = entrante[0] #Sale VB Saliente y entra la VB Entrante
    while i < len(nueva_matriz[0]):
        nueva_matriz[saliente[1]][i] = nueva_matriz[saliente[1]][i]/pivote
        i += 1
    nueva_matriz[saliente[1]][entrante[1]] = 1
    return nueva_matriz

def columna_seleccionada(matriz, entrante):
    """Retorna la columna pivote"""
    fila = 1
    columna = []
    while (fila < len(matriz)):
        columna += [-matriz[fila][entrante[1]]]
        fila += 1
    return columna

def iteracion(nueva_matriz, columna_iterada, fila_iterada, pos_columna_iterada, pos_pivote):
    num_fila = 0
    for fila in nueva_matriz:
        if num_fila != pos_pivote and num_fila != 0 :
            indice = 0
            for num in fila[1:]:

                if indice != pos_columna_iterada-1 and indice <= len(nueva_matriz) + 1:
                    nueva_matriz[num_fila][indice+1] = float(num+columna_iterada[num_fila-1]*fila_iterada[indice+1])
                indice += 1
                
               
        num_fila += 1     
    return nueva_matriz

def verificar_optimalidad(funcion_objetivo):
    for n in funcion_objetivo[1:]:
        if n < 0:
            return False
    return True

def encontrar_basicas(matriz):
    var_basicas = []
    for fila in matriz[2:]:
        var_basicas.append(fila[0])
    return var_basicas

def principal(args):
    
    diccionario_datos = {}
    matriz = []
    if len(args) == 3 and args[1] == "-h":
        diccionario_datos = leer_archivo(args[2])
        matriz = definir_ecuaciones(diccionario_datos)
        i = 0
        while(not(verificar_optimalidad(matriz[1]))):
            print("\n\nIteracion " + str(i))
            for fila in matriz:
                print(fila)
            FEV = encontrar_FEV(matriz)
            entrante = encontrar_entrante(matriz)
            saliente = encontrar_saliente(matriz, entrante)
            pivote = matriz[saliente[1]][entrante[1]]
            columna_iterada = columna_seleccionada(matriz, entrante)
            matriz = llenar_columna(matriz, entrante)
            matriz = llenar_fila(pivote, entrante, saliente, matriz)
            fila_iterada = matriz[saliente[1]]
            matriz = iteracion(matriz, columna_iterada, fila_iterada, entrante[1], saliente[1])
            print("FEV: " + str(FEV[1]))
            print("U: " + str(FEV[0]))
            print("Variable básica entrante: " + entrante[0])
            print("Variable básica saliente: " + saliente[0])
            print("Pivote: " + str(pivote))
            i += 1
        
        print("\n\nIteracion " + str(i))
        for fila in matriz:
            print(fila)
        FEV = encontrar_FEV(matriz)
        print("FEV: " + str(FEV[1]))
        print("U: " + str(FEV[0]))
        if FEV[2] != []:
            if FEV[2][0]:
                entrante = (FEV[2][1], FEV[2][2])
                saliente = encontrar_saliente(matriz, entrante)
                pivote = matriz[saliente[1]][entrante[1]]
                columna_iterada = columna_seleccionada(matriz, entrante)
                matriz = llenar_columna(matriz, entrante)
                matriz = llenar_fila(pivote, entrante, saliente, matriz)
                fila_iterada = matriz[saliente[1]]
                matriz = iteracion(matriz, columna_iterada, fila_iterada, entrante[1], saliente[1])
                print("Variable básica entrante: " + entrante[0])
                print("Variable básica saliente: " + saliente[0])
                print("Pivote: " + str(pivote))
                print("\n\nIteracion extra")
                for fila in matriz:
                    print(fila)
                FEV = encontrar_FEV(matriz)
                print("FEV: " + str(FEV[1]))
                print("U: " + str(FEV[0]))
    elif len(args) == 2:

        if args[1] == "-h":
            print("\nAca escribiremos el código de ayuda\n")
        
        else:
            diccionario_datos = leer_archivo(args[1])
            matriz = definir_ecuaciones(diccionario_datos)
            i = 0
            while(not(verificar_optimalidad(matriz[1]))):
                print("\n\nIteracion " + str(i))
                for fila in matriz:
                    print(fila)
                FEV = encontrar_FEV(matriz)
                entrante = encontrar_entrante(matriz)
                saliente = encontrar_saliente(matriz, entrante)
                pivote = matriz[saliente[1]][entrante[1]]
                columna_iterada = columna_seleccionada(matriz, entrante)
                matriz = llenar_columna(matriz, entrante)
                matriz = llenar_fila(pivote, entrante, saliente, matriz)
                fila_iterada = matriz[saliente[1]]
                matriz = iteracion(matriz, columna_iterada, fila_iterada, entrante[1], saliente[1])
                print("FEV: " + str(FEV[1]))
                print("U: " + str(FEV[0]))
                print("Variable básica entrante: " + entrante[0])
                print("Variable básica saliente: " + saliente[0])
                print("Pivote: " + str(pivote))
                i += 1
            
            print("\n\nIteracion " + str(i))
            for fila in matriz:
                print(fila)
            FEV = encontrar_FEV(matriz)
            print("FEV: " + str(FEV[1]))
            print("U: " + str(FEV[0]))
            if FEV[2] != []:
                if FEV[2][0]:
                    entrante = (FEV[2][1], FEV[2][2])
                    saliente = encontrar_saliente(matriz, entrante)
                    pivote = matriz[saliente[1]][entrante[1]]
                    columna_iterada = columna_seleccionada(matriz, entrante)
                    matriz = llenar_columna(matriz, entrante)
                    matriz = llenar_fila(pivote, entrante, saliente, matriz)
                    fila_iterada = matriz[saliente[1]]
                    matriz = iteracion(matriz, columna_iterada, fila_iterada, entrante[1], saliente[1])
                    print("Variable básica entrante: " + entrante[0])
                    print("Variable básica saliente: " + saliente[0])
                    print("Pivote: " + str(pivote))
                    print("\n\nIteracion extra")
                    for fila in matriz:
                        print(fila)
                    FEV = encontrar_FEV(matriz)
                    print("FEV: " + str(FEV[1]))
                    print("U: " + str(FEV[0]))

    else:
        print("\nIngrese [-h] para recibir ayuda de utilización del programa\n")
        return


principal(sys.argv)