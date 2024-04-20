#ifndef HASH_H
#define HASH_H

#include <iostream>
#include <cmath>
#include <random>
#include <chrono>

typedef std::size_t HASH_INDEX_T;

// my additions is this okay?
#include <string>
using namespace std;

struct MyStringHash {
    HASH_INDEX_T rValues[5] { 983132572, 1468777056, 552714139, 984953261, 261934300 };
    MyStringHash(bool debug = true)
    {
        if(false == debug){
            generateRValues();
        }
    }
    // hash function entry point (i.e. this is h(k))
    HASH_INDEX_T operator()(const std::string& k) const
    {
        // intialize helpful variables
        size_t k_length = k.length();
        HASH_INDEX_T res_arr[5];
        for(int i = 0; i < 5; i++){
          res_arr[i] = 0;
        }

        int res_index = 4;

        while((k_length / 6) != 0){
          // start with end of the string and extract 6-char substrings
          string substring = k.substr(k_length - 6, 6);
          // convert the characters in the substring to int value using helper
          HASH_INDEX_T arr[6];
          for(int i = 0; i < substring.length(); i++){
            arr[i] = letterDigitToNumber(substring[i]);
          }
          unsigned long long res = 0;
          // convert to a result unsigned long long
          for(int i = 0; i < 6; i++){
            res *= 36;
            res += arr[i];
          }
          // add to final result array and adjust indices
          res_arr[res_index] = res;
          res_index -= 1;
          k_length = k_length - 6;
          cout << k_length << endl;
        }

        // remaining chars that don't fit into 6-char substring
        int k_leftover = k_length % 6;

        // initialize leftover array to zeros
        HASH_INDEX_T arr[6];
        for(int i = 0; i < 6; i++){
          arr[i] = 0;
        }
        // extract leftover substring
        string substring = k.substr(0, k_length);
        // convert these letters to integer and add to array
        for(int i = 0; i < k_leftover; i++){
          arr[(6-k_leftover) + i] = letterDigitToNumber(substring[i]);
        }
        // convert arr to resulting unsigned long long
        unsigned long long res = 0;
        for(int i = 0; i < 6; i++){
          res *= 36;
          res += arr[i];
        }
        // add this to final result array
        res_arr[res_index] = res;
        
        // use values in res_arr to get hash
        HASH_INDEX_T hash_val = 0;
        for(int i = 0; i < 5; i++){
           hash_val += (res_arr[i] * rValues[i]);
        }
        // return final hash value
        return hash_val;
    }

    // A likely helper function is to convert a-z,0-9 to an integral value 0-35
    HASH_INDEX_T letterDigitToNumber(char letter) const
    {
        // Add code here or delete this helper function if you do not want it

        // convert to lower case if an upper case letter
        if (int(letter) >= 65 && int(letter) <= 90){
          letter = letter + 32;
        }
        // if in range a - z
        if(letter >= 'a' && letter <= 'z'){
          return letter - 'a';
        }
        // if one of the digits 0-9
        else {
         return letter - '0' + 26;
        }

    }

    // Code to generate the random R values
    void generateRValues()
    {
        // obtain a seed from the system clock:
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::mt19937 generator (seed);  // mt19937 is a standard random number generator

        // Simply call generator() [it has an operator()] to get another random number
        for(int i{ 0 }; i < 5; ++i)
        {
            rValues[i] = generator();
        }
    }
};

#endif
