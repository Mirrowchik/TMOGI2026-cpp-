#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <set>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;

struct Route {
    std::string whence;
    std::string whithere;
    double excess;
    double distance;
    std::string level;
};

struct routesimple {
    std::string name;
    double value;
};

std::vector<std::string> getUniqueNames(const std::vector<Route>& routes) {
    std::set<std::string> uniqueNames;
    for (const auto& route : routes) {
        uniqueNames.insert(route.whence);
        uniqueNames.insert(route.whithere);
    }
    return std::vector<std::string>(uniqueNames.begin(), uniqueNames.end());
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cerr << "Error: " << argv[0] << " <path.txt> <path.txt> <path.txt>" << std::endl;
        return 1;
    }

    std::ofstream resultFile("result.txt");
    if (!resultFile.is_open()) {
        std::cerr << "Error: Cannot create result.txt" << std::endl;
        return 1;
    }

    // Перенаправление std::cout в файл
    std::streambuf* coutBuffer = std::cout.rdbuf();
    std::cout.rdbuf(resultFile.rdbuf());


    // Чтение первого файла
    std::ifstream file1(argv[1]);
    if (!file1.is_open()) {
        std::cerr << "Error open: " << argv[1] << std::endl;
        return 1;
    }

    std::vector<Route> r1;
    std::string line;
    while (std::getline(file1, line)) {
        if (!line.empty()) {
            Route r;
            std::istringstream iss(line);
            std::getline(iss, r.whence, '\t');
            std::getline(iss, r.whithere, '\t');
            iss >> r.excess >> r.distance;
            std::getline(iss, r.level);
            r1.push_back(r);
        }
    }
    file1.close();
    VectorXd vExcess(r1.size());
    VectorXd vDistance(r1.size());
    vExcess.setZero();
    vDistance.setZero();
    for(int i = 0;  i<r1.size(); i++){vExcess(i)=r1[i].excess;vDistance(i)=r1[i].distance;}
    
    

    // Чтение второго файла
    std::ifstream file2(argv[2]);
    if (!file2.is_open()) {
        std::cerr << "Error open: " << argv[2] << std::endl;
        return 1;
    }

    std::vector<routesimple> r2;
    
    while (std::getline(file2, line)) {
        if (!line.empty()) {
            routesimple r;
            std::istringstream iss(line);
            std::getline(iss, r.name, '\t');  
            iss >> r.value;                   
            r2.push_back(r);
        }
    }
    file2.close();

    VectorXd vHigh(r2.size());
    vHigh.setZero();
    for(int i = 0;  i<r2.size(); i++){vHigh(i)+=r2[i].value; }
    
    // Получение уникальных имен и вычисление dof
    std::vector<std::string> uniqueNames = getUniqueNames(r1);
    const int n = static_cast<int>(r1.size());
    const int dof = static_cast<int>(r1.size() - uniqueNames.size() + r2.size());

    // Создание матрицы Eigen
    MatrixXd matrixB( dof,n);
    MatrixXd matrixF( dof,n);
    MatrixXd matrixHigh( dof, r2.size());
    // MatrixXd matrix(n + r2.size(), dof);
    matrixB.setZero();
    matrixF.setZero(); // Инициализация нулями
    matrixHigh.setZero();
    // matrix.setZero();

    // Чтение третьего файла в матрицу Eigen
    std::ifstream file3(argv[3]);
    if (!file3.is_open()) {
        std::cerr << "Error open: " << argv[3] << std::endl;
        return 1;
    }

    for (int j = 0; j < dof; ++j) {
        for (int i = 0; i < n + r2.size(); ++i) {
            double value;
            if (!(file3 >> value)) {
                std::cerr << "Error read matrixB " << n << "x" << dof << std::endl;
                return 1;
            }
            // Проверка допустимых значений
            if (value != -1.0f && value != 0.0f && value != 1.0f) {
                std::cerr << "invalid values " << value
                          << " Only allowed -1, 0, 1." << std::endl;
                return 1;
            }
            if(i<n){matrixB(j, i) = value;}
            else{matrixHigh( j,i-n) = value;}
        }
    }
    file3.close();

    std::ifstream file4(argv[4]);
    if (!file4.is_open()) {
        std::cerr << "Error open: " << argv[4] << std::endl;
        return 1;
    }

    for (int j = 0; j < dof; ++j) {
        for (int i = 0; i < n; ++i) {
            double value;
            if (!(file4 >> value)) {
                std::cerr << "Error read matrixB " << n << "x" << dof << std::endl;
                return 1;
            }
            // Проверка допустимых значений
            if (value != -1.0f && value != 0.0f && value != 1.0f) {
                std::cerr << "invalid values " << value
                          << " Only allowed -1, 0, 1." << std::endl;
                return 1;
            }
            matrixF(j, i) = value;
            
        }
    }
    file4.close();

    VectorXd vW(dof);
    vW.setZero();
    for(int i  = 0;i<dof;++i)
    {
        double sum = 0;
        for(int j=0;j<n;++j)
        {
            sum+= vExcess(j) * matrixB(i,j);
        }
        for(int j=0;j<r2.size();++j)
        {
            sum+= vHigh(j) * matrixHigh(i,j);
        }
        vW(i) = sum;
    }

    VectorXd vWd(dof);
    vWd.setZero();
    for(int i  = 0;i<dof;++i)
    {
        double sum = 0;
        for(int j=0;j<n;++j)
        {
            sum+= vDistance(j) * abs(matrixB(i,j));
        }
        if(r1[i].level == "\tIV"){vWd(i) =  20 * sqrt(sum);continue;}
        if(r1[i].level == "\tIII"){vWd(i) =  10 * sqrt(sum);continue;}
        if(r1[i].level == "\tII"){vWd(i) =  5 * sqrt(sum);continue;}
        
        
    }
    
    MatrixXd Q(n,n);
    Q.setZero();
    for(int i = 0;i<n;++i){Q(i,i)= 0.02  * sqrt(vDistance(i));}

    MatrixXd R(dof,dof);
    R.setZero();
    R = matrixB * Q * matrixB.transpose(); 

    VectorXd K(dof);
    K.setZero();
    K = -1 * R.inverse() * vW; 

    VectorXd V(n);
    V.setZero();
    V = (Q * matrixB.transpose()) * K;

    VectorXd BV(dof);
    BV.setZero();
    BV = matrixB * V;

    VectorXd vHighI(n);
    vHighI.setZero();
    vHighI=vExcess+V;

     VectorXd vHigha(dof);
    vHigha.setZero();
    for(int j = 0; j<n;j++){vHigha(0)+=vHighI(j)*matrixF(0,j);}vHigha(0)+=vHigh(0);
    for(int i = 1;i< dof;i++){for(int j = 0; j<n;j++){vHigha(i)+=vHighI(j)*matrixF(i,j);}vHigha(i)+=vHigh(1);}

    double u;
    {double sum = 0;for(int i = 0; i<n ; i++){sum += (1/Q(i,i)) * V(i) * V(i);}  u = sqrt(sum/dof);}

    MatrixXd Qv(n,n);
    Qv.setZero();
    Qv = -1 * Q * matrixB.transpose() * R.inverse() * matrixB * Q;
    
    MatrixXd Qy(n,n);
    Qy.setZero();
    Qy = Q + Qv; 

    VectorXd vmh(n);
    vmh.setZero();
    for(int i = 0;i< n;i++){vmh(i)=u*sqrt(Qy(i,i));}

    MatrixXd QH(dof,dof);
    QH.setZero();
    QH = matrixF * Qy * matrixF.transpose();

    VectorXd vmH(dof);
    vmH.setZero();
    for(int i = 0;i< dof;i++){vmH(i)=u*sqrt(QH(i,i));}

    // Вывод результатов
    std::cout << "dof = " << dof << std::endl;
    if (n > 0 && dof > 0) {
        std::cout << "matrixB " << std::endl << matrixB << std::endl;
        std::cout << "matrixH " << std::endl << matrixHigh << std::endl;
        std::cout << std::fixed << std::setprecision(7);
        std::cout << "VH " << std::endl << vHigh << std::endl;
        std::cout << "Vd " << std::endl << vDistance << std::endl;
        std::cout << "Ve " << std::endl << vExcess << std::endl;
        std::cout << "VW " << std::endl << vW << std::endl;
        std::cout << "VWd " << std::endl << vWd/1000 << std::endl;
        std::cout << "Q " << std::endl << Q.format(Eigen::IOFormat(0, Eigen::DontAlignCols, "\t│\t",/*между элементами*/"\n",/*между строками*/"│\t",/*начало строки*/"\t│",/*конец строки*/"", "")) << std::endl;
        std::cout << "R " << std::endl << R.format(Eigen::IOFormat(0, Eigen::DontAlignCols, "\t│\t",/*между элементами*/"\n",/*между строками*/"│\t",/*начало строки*/"\t│",/*конец строки*/"", "")) << std::endl;
        std::cout << "K " << std::endl << K << std::endl;
        std::cout << "V " << std::endl << V << std::endl;
        std::cout << "BV " << std::endl << BV << std::endl;
        std::cout << "vHighI " << std::endl << vHighI << std::endl;
        std::cout << "u " << std::endl << u << std::endl;
        std::cout << "Qv " << std::endl << Qv.format(Eigen::IOFormat(0, Eigen::DontAlignCols, "\t│\t",/*между элементами*/"\n",/*между строками*/"│\t",/*начало строки*/"\t│",/*конец строки*/"", "")) << std::endl;
        std::cout << "Qy " << std::endl << Qy.format(Eigen::IOFormat(0, Eigen::DontAlignCols, "\t│\t",/*между элементами*/"\n",/*между строками*/"│\t",/*начало строки*/"\t│",/*конец строки*/"", "")) << std::endl;
        std::cout << "vmh " << std::endl << vmh << std::endl;
        // std::cout << std::fixed << std::setprecision(0);
        std::cout << "matrixF " << std::endl << matrixF << std::endl;
        // std::cout << std::fixed << std::setprecision(7);
        std::cout << "QH " << std::endl << QH.format(Eigen::IOFormat(0, Eigen::DontAlignCols, "\t│\t",/*между элементами*/"\n",/*между строками*/"│\t",/*начало строки*/"\t│",/*конец строки*/"", "")) << std::endl;
        std::cout << "vmH " << std::endl << vmH << std::endl;
        std::cout << "vHigha " << std::endl << vHigha << std::endl;
    }

    // Восстановление стандартного вывода в консоль
    std::cout.rdbuf(coutBuffer);
    std::cout << "Results saved to result.txt" << std::endl;

    // Закрытие файла
    resultFile.close();

    // Дополнительные возможности Eigen (пример)
    // std::cout << "Matrix shape: " << matrixB.rows() << " x " << matrixB.cols() << std::endl;
    // std::cout << "Matrix norm: " << matrixB.norm() << std::endl;

    return 0;
}