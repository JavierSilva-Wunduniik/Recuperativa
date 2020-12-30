#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <mpi.h>
#include <algorithm>
#include <list>
#include <string>
#include <map>
#include <numeric>
#include <math.h>

std::map<double, double> smi;
std::map<double, double> smiaux;
std::map <std::string, double> dollars;
std::map <std::string, double> dollarsaux;
std::map<double, double> baccano;
std::list<double> lista;

auto convertYear(std::map<double, double> y, std::map<std::string, double> d) 
{
    std::map<double, double> baccano;
    auto anos = y.begin(); // Inicializa variable e iteracion       
    while (anos != y.end()) 
    { 
        // Mientras "anos" no haya terminado este ciclo se mantiene en ejecución
        lista = {}; // Reincia lista por cada año
        auto dollars = d.begin();
        while (dollars != d.end()) 
        {            
            if (anos->first == std::stod(dollars->first.substr(0, 4))) 
            { 
                //Comparacion entre año del SMI y el año de los valores diarios del dollar                
                lista.push_back(dollars->second); // Agrega los valores de dollar a la lista
            }
            ++dollars;
        }
        double suma = std::accumulate(std::begin(lista), std::end(lista), 0.0); // Suma todos los valores de la lista de los valores de dollar
        double aux = suma / lista.size(); //promedio
        baccano.insert(std::pair<double, double>(anos->first, anos->second/aux)); //salario minimo convertido a dollar
        ++anos;
    }
    return baccano;
}
auto linearRegression(std::map<double, double> m) 
{
    double n = m.size(); // m (mapa)
    double sx = 0;
    double sxx = 0;
    double sy = 0;
    double syy = 0;
    double sxy = 0;
    auto iter = m.begin(); // Inicializando la iteracion del mapa
    while (iter != m.end()) 
    {
        sx += iter->first;                 //sumatoria años salario minimo
        sxx += pow(iter->first, 2);        // x**2 - SUMA
        sy += iter->second;                //sumatoria salario minimo en dollar
        syy += pow(iter->second, 2);       // y**2 - SUMA
        sxy += iter->first * iter->second; // Multiplica año por valor corresp del dollar
        ++iter;
    }
    double b = (n * sxy - sx * sy) / (n * sxx - sx * sx); // Pendiente
    double a = (1 / n) * sy - b * (1 / n) * sx; // Coeficiente
    std::cout << std::endl;
    std::cout <<"====== Resultado ======"<< std::endl;
    std::cout << std::endl;
    std::cout << "y = " << b << "x " << a;
}

// Esta función revisa el archivo smi y obtiene sus valores para ingresarlos a otra variable con la que si se puede trabajar
auto csv_smi(char** argv) 
{
    std::map<double, double> map_tmp;
    std::ifstream csv;
    csv.open(argv[1]);
    std::string line;
    getline(csv, line); // Se mite la primera linea de descripcion
    while (csv.peek() != EOF) 
    {        
        getline(csv, line);
        line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());     // Se eliminan los saltos de linea
        size_t position = line.find(";") + 1;                                    // Guarda las posiciones de ;
        std::string key_k = line.substr(0, line.find(";"));                      // Esta variable toma lo que este antes de ;
        std::string key_v = line.substr(position);                               // En cambio esta toma lo que está después de ;
        key_k.erase(std::remove(key_k.begin(), key_k.end(), '\"'), key_k.end()); // Se eliminan las comillas
        key_v.erase(std::remove(key_v.begin(), key_v.end(), '\"'), key_v.end()); // Se eliminan las comillas
        double key_kd = std::stod(key_k);                                        // Se transforman las variables a double
        double key_vd = std::stod(key_v);
        map_tmp.insert(std::pair<double, double>(key_kd, key_vd));              // Se ingresan solo los valores que SI se quieren del archivo de forma que el programa pueda procesarlos
    }
    csv.close();
    return map_tmp;
}

// Esta función ingresa al archivo dollars y obtiene los valores del mismo para ingresarlos a una variable de forma que si se puedan trabajar.
auto csv_dollars(char** argv) 
{
    std::map<std::string, double> map_tmp;
    std::ifstream csv;
    csv.open(argv[2]);
    std::string line;
    getline(csv, line); // Se omite la primera linea de descripcion
    while (csv.peek() != EOF) 
    { 
        //Revisa si ha terminado de leer el archivo
        getline(csv, line);                                                      // Guarda linea del archivo en variable line
        line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());     // Elimina los saltos de linea
        size_t position = line.find(";") + 1;                                    // Guarda las posiciones de ;
        std::string key_k = line.substr(0, line.find(";"));                      // Variable toma lo que este antes de ;
        std::string key_v = line.substr(position);                               // En cambio esta toma lo que está después de ;
        key_k.erase(std::remove(key_k.begin(), key_k.end(), '\"'), key_k.end()); // Se eliminan las comillas
        key_v.erase(std::remove(key_v.begin(), key_v.end(), '\"'), key_v.end()); // Se eliminan las comillas
        double key_vd = std::stod(key_v);                                        // Convierte string valor dollar en doble
        map_tmp.insert(std::pair<std::string, double>(key_k, key_vd));   // Se ingresan solo los valores que SI se quieren del archivo de forma que el programa pueda procesarlos
    }
    csv.close();
    return map_tmp;    
}

//Este es el "main" 
int main(int argc, char** argv) 
{

    int nodo, comuni;
    MPI_Init(&argc,&argv);
    MPI_Status status;
    MPI_Comm_rank( MPI_COMM_WORLD , &nodo);
    MPI_Comm_size(MPI_COMM_WORLD, &comuni);
    if (nodo==0)
    {
        smi = csv_smi(argv);
        MPI_Send(&smi,0,MPI_DOUBLE ,2,0,MPI_COMM_WORLD);
    }
        if (nodo==1)
    {
        dollars = csv_dollars(argv);
        MPI_Send(&dollars,0,MPI_DOUBLE ,2,0,MPI_COMM_WORLD);
    }
    if (nodo==2)
    {
        MPI_Recv(&smiaux,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&status);
        MPI_Recv(&dollarsaux,1,MPI_DOUBLE,1,0,MPI_COMM_WORLD,&status);
        baccano = convertYear(smiaux, dollarsaux);
        linearRegression(baccano);
    }
    MPI_Finalize;
    /*
    smi = csv_smi(argv);
    dollars = csv_dollars(argv);
    baccano = convertYear(smi, dollars);
    linearRegression(baccano);
    */
    
    std::cout << std::endl;
    std::cout <<"====== Integrante ======"<< std::endl;
    std::cout << std::endl;
    std::cout << "Javier Silva G"<< std::endl;
    std::cout << std::endl;
    return 0;
}