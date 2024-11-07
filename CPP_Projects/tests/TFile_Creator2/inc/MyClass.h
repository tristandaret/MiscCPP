#include <TH1F.h>
#include <TObject.h>

class MyClass : public TObject {
public:
    ClassDef(MyClass, 1);  // ROOT versioning for I/O

	MyClass();
    
    virtual ~MyClass();

    void FillHistogram(double value);

	TH1F *GetHistogram() const { return ph1f; }	

	private:
    TH1F *ph1f;  // Pointer to a TH1F histogram
};
