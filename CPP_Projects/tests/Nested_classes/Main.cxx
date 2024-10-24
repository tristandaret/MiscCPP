#include "Outer.h"

int main() {
    Outer outer;
    Outer::Inner inner(outer);  // Pass the Outer instance to Inner's constructor
    inner.accessOuterPrivateFields();  // Now, we don't need to pass the outer object

    return 0;
}