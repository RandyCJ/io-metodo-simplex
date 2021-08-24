import sys

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
    print(nueva_matriz)

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
    crear_matriz([diccionario_datos["fun_ob"]] + diccionario_datos["rest"], var_basicas, diccionario_datos["num_var"] + diccionario_datos["num_rest"])

def principal(args):
    
    diccionario_datos = {}

    if len(args) == 3 and args[1] == "-h":
        print ("\nExiste el argumento de ayuda y el argumento del archivo\n")
        diccionario_datos = leer_archivo(args[2])
        definir_ecuaciones(diccionario_datos)
    
    elif len(args) == 2:

        if args[1] == "-h":
            print("\nAca escribiremos el código de ayuda\n")
        
        else:
            print("\nSolo el argumento del archivo\n")
            diccionario_datos = leer_archivo(args[1])
            definir_ecuaciones(diccionario_datos)

    else:
        print("\nIngrese [-h] para recibir ayuda de utilización del programa\n")
        return


principal(sys.argv)