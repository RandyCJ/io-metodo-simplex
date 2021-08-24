import sys

def leer_archivo(nombre_archivo):
    
    contador = 0
    lista_datos = []
    lista_tmp = []
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

def principal(args):
    
    if len(args) == 3 and args[1] == "-h":
        print ("\nExiste el argumento de ayuda y el argumento del archivo\n")
        leer_archivo(args[2])
    
    elif len(args) == 2:

        if args[1] == "-h":
            print("\nAca escribiremos el código de ayuda\n")
        
        else:
            print("\nSolo el argumento del archivo\n")
            leer_archivo(args[1])

    else:
        print("\nIngrese [-h] para recibir ayuda de utilización del programa\n")


principal(sys.argv)