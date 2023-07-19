# SGA
This is a C++ geometric algebra library intended to show how simple it is to create a general GA library without being that complicated.  The one header file is less than 600 lines of code, and a decent amount of it is just boilerplate to make it slightly easier to use.  The code is fully general, allowing for any dimension (technically up to 64) and any metric.  The code also only keeps track of multivector coefficients when needed by using templates.  Furthermore, very little template metaprogramming is needed.

This library is mainly intended for educational purposes.  While it should be decently efficient, know that it hasn't really been optimized to squeeze out every last drop of performance.

If you do wish to use this library in your code, it's just a header-only library, so feel free to put the header into your project.  The code is in the public domain so you can do what you wish with it.  The code uses many C++20 features, so you'll need a modern compiler to get it to work.
