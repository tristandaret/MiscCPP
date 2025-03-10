#include <iostream>
#include <TMatrixD.h>
#include <TRandom3.h>
#include <iomanip>

int main()
{
   // Initialize random generator
   TRandom3 randGen(42); // Seed for reproducibility

   const double scaleFactor = 1e15;
   // Create a 5x5 symmetric matrix
   TMatrixD mat(5, 5);
   for (int i = 0; i < 5; ++i) {
      for (int j = i; j < 5; ++j) {
         double value = randGen.Uniform(-10, 10) * scaleFactor; // Random values between -10 and 10
         std::cout << std::setprecision(15) << value << std::endl;
         mat(i, j) = value;
         mat(j, i) = value; // Ensure symmetry
      }
   }

   // Compute the inverse
   TMatrixD matInv = mat;
   matInv.Invert();

   // Multiply matrix by its inverse
   TMatrixD identity = mat * matInv;

   // Print the matrices
   std::cout << "Matrix:";
   mat.Print();
   std::cout << "Inverse:";
   matInv.Print();
   std::cout << "Matrix * Inverse:";
   identity.Print();

   // Check if the result is a diagonal matrix
   bool isDiagonal = true;
   for (int i = 0; i < 5; ++i) {
      for (int j = 0; j < 5; ++j) {
         if (i != j && std::abs(identity(i, j)) > 1e-10) { // Allow small numerical errors
            isDiagonal = false;
            break;
         }
      }
      if (!isDiagonal)
         break;
   }

   // Print the result of the check
   if (isDiagonal) {
      std::cout << "The result is a diagonal matrix (identity matrix)." << std::endl;
   } else {
      std::cout << "The result is NOT a perfect diagonal matrix (numerical errors possible)." << std::endl;
   }

   return 0;
}
