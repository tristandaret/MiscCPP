#include "Outer.h"

// Outer constructor (default constructor)
Outer::Outer() = default;

// Inner class constructor implementation
Outer::Inner::Inner(Outer &o) : outer(o) {}

// Method to access Outer class private fields
void Outer::Inner::accessOuterPrivateFields()
{
   std::cout << "privateField1: " << outer.privateField1 << std::endl;
   std::cout << "privateField2: " << outer.privateField2 << std::endl;
   std::cout << "privateField3: " << outer.privateField3 << std::endl;
}
