#ifndef OUTER_H
#define OUTER_H

#include <iostream>
#include <string>

class Outer {
private:
   int privateField1 = 0;                 // Default value
   std::string privateField2 = "default"; // Default value
   double privateField3 = 3.14;           // Default value

   friend class Inner; // Inner is a friend of Outer and can access private fields

public:
   Outer(); // Declaration of default constructor

   class Inner {
   private:
      Outer &outer; // Reference to Outer class

   public:
      Inner(Outer &o); // Constructor declaration

      void accessOuterPrivateFields(); // Method declaration
   };
};

#endif // OUTER_H
