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
    var_holgura = [] 
    soluciones_multiples = False
    nom_archivo = ""
    U = 0
    FEV = []
    dual = False
    CONST_M = 0
    var_artificiales = []
    es_max = True

    def __init__(self, matriz) -> None:
        self.matriz = matriz
        if matriz[1][0] == "-U":
            self.es_max = False
    
    def definir_artificiales(self, var_artificiales, M):
        """ Define las variables artificiales en la matriz con A1, A2...
            E: las variables artificiales en orden de X4, X6.., la constante M
            S: N/A
        """
        artificiales_tmp = []
        for n in var_artificiales:
            artificiales_tmp.append("X" + str(n))

        i = 1
        num_artificial = 1
        num_normal = 1
        # Primero cambia las variables en la primera fila
        while i < len(self.matriz[0])-1:
            if self.matriz[0][i] in artificiales_tmp:
                self.matriz[0][i] = "A" + str(num_artificial)
                self.var_artificiales.append("A" + str(num_artificial))
                num_artificial += 1
            else:
                self.matriz[0][i] = "X" + str(num_normal)
                num_normal += 1
            i += 1

        #Ahora las cambia de la columna de variables basicas
        i = 2
        num_artificial = 1
        while i < len(self.matriz):
            if self.matriz[i][0] in artificiales_tmp:
                self.matriz[i][0] = "A" + str(num_artificial)
                num_artificial += 1
            i += 1
        self.CONST_M = M

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
            raise Exception("No acotada")

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
        if self.es_max:
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
                matriz_redondeada[fila][columna] = self.matriz[fila][columna]
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
        if self.es_max:
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

    def verificar_artificiales(self):
        """ Verifica en la solucion óptima que las variables artificiales no sean positivas
            que significaría que la solución es no factible
            E: N/A
            S: Retorna 0 si hay solucion factible, caso contrario retorna el X artificial que seria un string
        """
        fila = 2

        while fila < len(self.matriz)-1:
            if self.matriz[fila][-1] > 0 and self.matriz[fila][0] in self.var_artificiales:
                return self.matriz[fila][0]
            fila += 1
        return 0