#include <iostream>
using namespace std;

int* doubler(int a[], int size)
{
    int b[size];
    for (int i=0; i<size+1; i++) b[i]=2*a[i];
    int* pointer = b;
    return pointer;
}

int main()
{
    int a[5]={0};
//    for (int i=0;i<5;i++) printf("%d ",a[i]); printf("\n");
//    printf("a[4] = %d\n",a[4]);
    int* b;
    b = doubler(a,5);
    for (int i=0;i<5;i++) printf("%d ",b[i]);
    printf("\n");
    return 0;
}
