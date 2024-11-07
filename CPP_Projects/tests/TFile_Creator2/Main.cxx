#include <iostream>
#include <TFile.h>
#include "MyClass.h"

void SaveMyClassToFile() {
    // Create an instance of the custom class
    MyClass* myObject = new MyClass();
    myObject->FillHistogram(2.3);  // Example of filling the histogram

    // Open a ROOT file
    TFile file("myFile.root", "RECREATE");

    // Write the object to the file
    file.WriteObject(myObject, "myObject");

    // Close the file
    file.Close();

    delete myObject;  // Clean up
}

void LoadMyClassFromFile() {
    TFile file("myFile.root", "READ");
    MyClass* myObject;
    file.GetObject("myObject", myObject);

    if (myObject && myObject->GetHistogram()) {
        std::cout << "mean: " << myObject->GetHistogram()->GetMean() << std::endl;
    }
    
    delete myObject;  // Clean up
}

int main() {
	SaveMyClassToFile();
	LoadMyClassFromFile();
	return 0;
}