#include <iostream>
#include <cmath>
#include <vector>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>
#include <iostream>

using namespace std;

// Класс расчётной точки
class CalcNode
{
// Класс сетки будет friend-ом точки
friend class CalcMesh;

protected:
    // Координаты
    double x;
    double y;
    double z;
    // Некая величина, в попугаях
    double smth;
    // Скорость
    double vx;
    double vy;
    double vz;

public:
    // Конструктор по умолчанию
    CalcNode() : x(0.0), y(0.0), z(0.0), vx(0), vy(0), vz(0)
    {
    }

    // Конструктор с указанием всех параметров
    CalcNode(double x, double y, double z, double vx, double vy, double vz) 
            : x(x), y(y), z(z), vx(vx), vy(vy), vz(vz)
    {
    }
};

// Класс расчётной сетки
class CalcMesh
{
protected:
    // 2D-сетка из расчётных точек
    vector<vector<CalcNode>> points;
public:
    // Конструктор сетки size x size точек с шагом h по пространству
    CalcMesh(unsigned int size, double h) {
        points.resize(size);
        for(unsigned int i = 0; i < size; i++) {
            points[i].resize(size);
            if (i < size - 1) points[i+1].resize(size);
            for(unsigned int j = 0; j < size; j++) {
            double side = size * h;

                // Начальные координаты зададим равномерно в плоскости OXY
                double pointX = i * h;
                double pointY = j * h;
                double pointZ = 0;

                double x = pointX - side / 2;
                double y = pointY - side / 2;

                if (x*x + y*y < 0.16) {
                    pointZ = 5 * (0.16 - x*x - y*y);
                }

                // pointZ = 0.5 * sin(6.28 * pointX * 3) * sin (6.28 * pointY * 3);
                
                points[i][j] = CalcNode(pointX, pointY, pointZ, 0, 0, 0);
            }
        }
    }

    // Метод отвечает за выполнение для всей сетки шага по времени величиной tau
    void doTimeStep(double tau, double h, int size) {
        // По сути метод просто двигает все точки
        for(unsigned int i = 0; i < size; i++) {
            for(unsigned int j = 0; j < size; j++) {

                double c = 0.5;

                double z0 = points[i][j].z;
                double v0 = points[i][j].vz;

                auto compute_acceleration = [&](int idx_i, int idx_j) -> double {
                    // Проверка границ (для краёв сетки)
                    if (idx_i == 0 || idx_i >= size - 1 ||
                        idx_j == 0 || idx_j >= size - 1) {
                        return 0.0;  // Граничные условия: закреплённая граница
                    }
                    
                    // Вычисление вторых производных (лапласиан)
                    double d2z_dx2 = (points[idx_i][idx_j+1].z - 2*points[idx_i][idx_j].z + points[idx_i][idx_j-1].z) / (h*h);
                    double d2z_dy2 = (points[idx_i+1][idx_j].z - 2*points[idx_i][idx_j].z + points[idx_i-1][idx_j].z) / (h*h);
                    
                    return c * c * (d2z_dx2 + d2z_dy2);
                };

                double k1_z = v0;
                double k1_v = compute_acceleration(i, j);
                
                // Сохраняем исходные значения для промежуточных вычислений
                double original_z = points[i][j].z;
                double original_v = points[i][j].vz;
                
                // k2: смещаем на полшага
                points[i][j].z = z0 + 0.5 * tau * k1_z;
                points[i][j].vz = v0 + 0.5 * tau * k1_v;
                double k2_z = points[i][j].vz;
                double k2_v = compute_acceleration(i, j);
                
                // k3: ещё один полшага
                points[i][j].z = z0 + 0.5 * tau * k2_z;
                points[i][j].vz = v0 + 0.5 * tau * k2_v;
                double k3_z = points[i][j].vz;
                double k3_v = compute_acceleration(i, j);
                
                // k4: полный шаг
                points[i][j].z = z0 + tau * k3_z;
                points[i][j].vz = v0 + tau * k3_v;
                double k4_z = points[i][j].vz;
                double k4_v = compute_acceleration(i, j);
                
                // Финальное обновление
                points[i][j].z = z0 + tau / 6.0 * (k1_z + 2*k2_z + 2*k3_z + k4_z);
                points[i][j].vz = v0 + tau / 6.0 * (k1_v + 2*k2_v + 2*k3_v + k4_v);
            }
        }
    }

    // Метод отвечает за запись текущего состояния сетки в снапшот в формате VTK
    void snapshot(unsigned int snap_number) {
        // Сетка в терминах VTK
        vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
        // Точки сетки в терминах VTK
        vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

        // Скалярное поле на точках сетки
        auto smth = vtkSmartPointer<vtkDoubleArray>::New();
        smth->SetName("smth");

        // Векторное поле на точках сетки
        auto vel = vtkSmartPointer<vtkDoubleArray>::New();
        vel->SetName("velocity");
        vel->SetNumberOfComponents(3);

        // Обходим все точки нашей расчётной сетки
        unsigned int number = (unsigned int)points.size();
        for(unsigned int i = 0; i < number; i++) {
            for(unsigned int j = 0; j < number; j++) {
                // Вставляем новую точку в сетку VTK-снапшота
                dumpPoints->InsertNextPoint(points[i][j].x, points[i][j].y, points[i][j].z);

                // Добавляем значение векторного поля в этой точке
                double _vel[3] = {points[i][j].vx, points[i][j].vy, points[i][j].vz};
                vel->InsertNextTuple(_vel);

                // И значение скалярного поля тоже
                smth->InsertNextValue(points[i][j].smth);
            }
        }

        // Задаём размеры VTK-сетки (в точках, по трём осям)
        structuredGrid->SetDimensions(number, number, 1);
        // Грузим точки в сетку
        structuredGrid->SetPoints(dumpPoints);

        // Присоединяем векторное и скалярное поля к точкам
        structuredGrid->GetPointData()->AddArray(vel);
        structuredGrid->GetPointData()->AddArray(smth);

        // Создаём снапшот в файле с заданным именем
        string fileName = "anim/membrane-step-" + std::to_string(snap_number) + ".vts";
        vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
        writer->SetFileName(fileName.c_str());
        writer->SetInputData(structuredGrid);
        writer->Write();
    }
};

int main()
{
    // Размер расчётной сетки, точек на сторону
    unsigned int size = 100;
    // Шаг точек по пространству
    double h = 0.01;
    // Шаг по времени
    double tau = 0.01;

    // Создаём сетку заданного размера
    CalcMesh mesh(size, h);

    // Пишем её начальное состояние в VTK
    mesh.snapshot(0);

    // Делаем шаги по времени, 
    // на каждом шаге считаем новое состояние и пишем его в VTK
    for(unsigned int step = 1; step < 250; step++) {
        mesh.doTimeStep(tau, h, size);
        mesh.snapshot(step);
    }

    return 0;
}
