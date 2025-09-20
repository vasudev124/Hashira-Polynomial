#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>

#include <gmp.h>

#include "json.hpp"

using json = nlohmann::json;
using namespace std;

void decodeY(mpz_t result, const string& value_str, int base) {
    if (mpz_set_str(result, value_str.c_str(), base) != 0) {
        throw runtime_error("Failed to decode Y value: '" + value_str + "' with base " + to_string(base));
    }
}

void findSecret(mpz_t secret_c, const vector<long>& x_coords, const vector<mpz_t>& y_coords) {
    mpz_set_ui(secret_c, 0);

    size_t k = x_coords.size();
    
    mpz_t term, numerator, denominator;
    mpz_init(term);
    mpz_init(numerator);
    mpz_init(denominator);

    for (size_t i = 0; i < k; ++i) {
        mpz_set_ui(numerator, 1);
        mpz_set_si(denominator, 1);

        for (size_t j = 0; j < k; ++j) {
            if (i == j) continue;
            mpz_mul_si(numerator, numerator, -x_coords[j]);
            mpz_mul_si(denominator, denominator, x_coords[i] - x_coords[j]);
        }
        
        mpz_mul(term, y_coords[i], numerator);
        mpz_divexact(term, term, denominator);

        mpz_add(secret_c, secret_c, term);
    }
    
    mpz_clear(term);
    mpz_clear(numerator);
    mpz_clear(denominator);
}

void solve_for_file(const string& filename) {
    ifstream file_stream(filename);
    if (!file_stream.is_open()) {
        throw runtime_error("Could not open file: " + filename);
    }
    json data = json::parse(file_stream);

    int k = data["keys"]["k"];
    vector<long> x_coords;
    vector<mpz_t> y_coords;
    
    y_coords.resize(k);
    for(int i = 0; i < k; ++i) {
        mpz_init(y_coords[i]);
    }

    int count = 0;
    for (auto& element : data.items()) {
        if (count >= k) break;
        if (element.key() == "keys") continue;

        x_coords.push_back(stol(element.key()));
        
        int base = stoi(element.value()["base"].get<string>());
        string value_str = element.value()["value"].get<string>();

        decodeY(y_coords[count], value_str, base);
        count++;
    }

    mpz_t secret_c;
    mpz_init(secret_c);
    findSecret(secret_c, x_coords, y_coords);

    char* secret_str = mpz_get_str(NULL, 10, secret_c);
    cout << secret_str << endl;
    
    free(secret_str);
    mpz_clear(secret_c);
    for (auto& y : y_coords) {
        mpz_clear(y);
    }
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Error: Please provide two JSON file paths as arguments." << endl;
        cerr << "Usage: ./solver testcase1.json testcase2.json" << endl;
        return 1;
    }

    try {
        solve_for_file(argv[1]);
        solve_for_file(argv[2]);
    } catch (const exception& e) {
        cerr << "An error occurred: " << e.what() << endl;
        return 1;
    }
    return 0;
}


