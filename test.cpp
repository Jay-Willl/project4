#include <random>
#include <chrono>

#include "matrix.hpp"

int randInt(int min, int max) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(min, max);
    return dis(gen);
}

float randFloat(float min, float max) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);
    return (float)dis(gen);
}

double randDouble(double min, double max) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);
    return dis(gen);
}

void testBasic() {
    Matrix<int> m1(3, 3, 1);
    int arr1[] = {0,1, 2, 3, 4, 5, 6, 7, 8, 9};
    m1.setElements(arr1);
    std::cout << m1 << std::endl;


    Matrix<int> m2(128, 128, 3);
    int len2 = 128 * 128 * 3;
    int *arr2 = new int[len2];
    for (int i = 0; i < len2; i++) {
        arr2[i] = randInt(0, 1000);
    }
    m2.setElements(arr2);
//    std::cout << m2 << std::endl;

    std::cout << m1.getElementValue(1, 1, 0) << std::endl;
//    std::cout << m1.getElementValue(1, 1, 1) << std::endl;
    std::cout << m2.getElementValue(100, 100, 0) << std::endl;
    std::cout << m2.getElementValue(100, 100, 1) << std::endl;
    std::cout << m2.getElementValue(100, 100, 2) << std::endl;
}

void testOperator() {
    Matrix<int> mi1(3, 2, 1);
    int arr1[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    mi1.setElements(arr1);

}

void testMul() {
    Matrix<int> mi1(3, 2, 1);
    int arr1[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    mi1.setElements(arr1);
    std::cout << mi1 << std::endl;
    std::cout << "----------------" << std::endl;
    Matrix<int> mi2(2, 3, 1);
    mi2.setElements(arr1);
    std::cout << mi2 << std::endl;
    std::cout << "----------------" << std::endl;
    Matrix<int> mi3 = mi2 * 3;
    std::cout << mi3 << std::endl;
    std::cout << "----------------" << std::endl;
    Matrix<int> mi4(3, 3, 1);
    multiply(mi1, mi2, mi4);
    std::cout << "----------------" << std::endl;

    Matrix<int> mi5;
    auto *arri1 = new int[mi1.total];
    for (int i = 0; i < mi1.total; ++i) {
        arri1[i] = randInt(0, 100);
    }
    Matrix<int> mi6;
    auto *arri2 = new int[mi2.total];
    for (int i = 0; i < mi2.total; ++i) {
        arri2[i] = randInt(0, 100);
    }
    Matrix<int> mi7;
    std::cout << ">>>>>>>>>>>>>>>>" << std::endl;
    auto total_i = 0;
    for (int i = 0; i < 1000; ++i) {
        auto start_time_i = std::chrono::high_resolution_clock::now();
        multiply(mi5, mi6, mi7);
        auto end_time_i = std::chrono::high_resolution_clock::now();
        auto duration_i = std::chrono::duration_cast<std::chrono::microseconds>(end_time_i - start_time_i).count();
        std::cout << "TEST: multiply between int round " << i << " takes " << duration_i << " ms" << std::endl;
        total_i += duration_i;
    }
    std::cout << "TEST: multiply between int matrices takes " << total_i / 1000 << " ms (avg)" << std::endl;
    std::cout << ">>>>>>>>>>>>>>>>" << std::endl;
    std::cout << "----------------" << std::endl;

    Matrix<float> mf1;
    auto *arrf1 = new float[mf1.total];
    for (int i = 0; i < mf1.total; ++i) {
        arrf1[i] = randFloat(0.0f, 100.f);
    }
    Matrix<float> mf2;
    auto *arrf2 = new float[mf2.total];
    for (int i = 0; i < mf2.total; ++i) {
        arrf2[i] = randFloat(0.0f, 100.f);
    }
    Matrix<float> mf3;
    std::cout << ">>>>>>>>>>>>>>>>" << std::endl;
    auto total_f = 0;
    for (int i = 0; i < 1000; ++i) {
        auto start_time_f = std::chrono::high_resolution_clock::now();
        multiply(mf1, mf2, mf3);
        auto end_time_f = std::chrono::high_resolution_clock::now();
        auto duration_f = std::chrono::duration_cast<std::chrono::microseconds>(end_time_f - start_time_f).count();
        std::cout << "TEST: multiply between float round " << i << " takes " << duration_f << " ms" << std::endl;
        total_f += duration_f;
    }
    std::cout << "TEST: multiply between float matrices takes " << total_f / 1000 << " ms (avg)" << std::endl;
    std::cout << ">>>>>>>>>>>>>>>>" << std::endl;
    std::cout << "----------------" << std::endl;

    Matrix<double> md1;
    auto *arrd1 = new double[md1.total];
    for (int i = 0; i < md1.total; ++i) {
        arrd1[i] = randDouble(0.0, 100.0);
    }
    Matrix<double> md2;
    auto *arrd2 = new double[md2.total];
    for (int i = 0; i < md2.total; ++i) {
        arrd2[i] = randDouble(0.0, 100.0);
    }
    Matrix<double> md3;
    std::cout << ">>>>>>>>>>>>>>>>" << std::endl;
    auto total_d = 0;
    for (int i = 0; i < 1000; ++i) {
        auto start_time_d = std::chrono::high_resolution_clock::now();
        multiply(md1, md2, md3);
        auto end_time_d = std::chrono::high_resolution_clock::now();
        auto duration_d = std::chrono::duration_cast<std::chrono::microseconds>(end_time_d - start_time_d).count();
        std::cout << "TEST: multiply between double round " << i << " takes " << duration_d << " ms" << std::endl;
        total_d += duration_d;
    }
    std::cout << "TEST: multiply between double matrices takes " << total_d / 1000 << " ms (avg)" << std::endl;
    std::cout << ">>>>>>>>>>>>>>>>" << std::endl;
    std::cout << "----------------" << std::endl;

}

int main() {
    std::cout << "just for test" << std::endl;
    testBasic();
    testMul();
}
