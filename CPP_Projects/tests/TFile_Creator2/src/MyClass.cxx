#include "MyClass.h"

ClassImp(MyClass);

// Default constructor
MyClass::MyClass(){};

// Destructor
MyClass::~MyClass()
{
   delete ph1f;
}

void MyClass::FillHistogram(double value)
{
   ph1f->Fill(value);
}
