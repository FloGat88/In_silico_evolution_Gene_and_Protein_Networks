//
//  Utility.hpp
//  EmbrionicDevelopment
//
//  Created by Florian Gartner on 22/12/16.
//  Copyright © 2016 Florian Gartner. All rights reserved.
//
// contains generally utile functions for Embryonic Development



#ifndef Utility_hpp
#define Utility_hpp

#include <stdio.h>
#include <iostream>
#include <vector>

using namespace std;

// define a Macro that serves as Tic Toc Token equivalently to that of Matlab
#define _TIC_ auto begin = std::chrono::high_resolution_clock::now()
#define _TOC_ auto end = std::chrono::high_resolution_clock::now(); std::cout << std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count()*1.0e-6 << "sec" << std::endl


template<class T>
void VecDel_element(vector<T>& vec, const T& val)          // VectorDelete: in vector vec find (the last) element that is identical to val and delete it, i.e. transfer the last vector element on its position and pop_back
{
    typename vector<T>::iterator it;
    it = vec.end()-1;
    while (*it != val && it>=vec.begin())                              // search vector from behind until we find an element that is identical to val
        --it;
    
    if(it<vec.begin()) {        // !!! Control mechanism: to be removed later !!!
        cout << "Element not found in vec by VecDel_element" << endl;
        return;
    }
    if(it!=vec.end()-1)                     // transfer last element to that position if it is not the last position anyways.
        *it = vec.back();
    vec.pop_back();                         // delete last element via pop_back
}



#endif /* Utility_hpp */
