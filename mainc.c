#include <stdio.h>
#include <gmp.h>
#include<stdlib.h>
#include<string.h>

typedef struct {
    mpz_t x;
    mpz_t y;
}Coordinates ;

void print_coordinates(const Coordinates *coord) {
    // Conversion de mpz_t en chaînes de caractères
    char *x_str = mpz_get_str(NULL, 10, coord->x);
    char *y_str = mpz_get_str(NULL, 10, coord->y);

    // Affichage des coordonnées
    printf("Coordonnees : x = %s, y = %s\n", x_str, y_str);

    // Libération de la mémoire allouée par mpz_get_str
    free(x_str);
    free(y_str);
}

Coordinates double_point(mpz_t x, mpz_t y ,mpz_t a ,mpz_t p){

    Coordinates result;

    mpz_inits(result.x , result.y, NULL);

    if (x == 0 && y ==0){
        
        mpz_set_ui(result.x, 0);
        mpz_set_ui(result.y, 0);
        return result;
    }

    mpz_t x3;
    mpz_init(x3);
    mpz_mul_si(x3, x,3);

    mpz_t two_mul_y;
    mpz_init(two_mul_y);
    mpz_mul_si(two_mul_y, y,2);

    mpz_t invert;
    mpz_init(invert);
    mpz_invert(invert, two_mul_y, p);


    mpz_t s; //we want s = ((3 * (x *x) + a ) * invert) % p;
    mpz_init(s);

    mpz_t x_square;
    mpz_init(x_square);
    mpz_mul(x_square, x, x);

    mpz_t three_x_square;
    mpz_init(three_x_square);
    mpz_mul_si(three_x_square,x_square,3);

    mpz_t the_up_of_the_fraction_aka_numerator;
    mpz_init(the_up_of_the_fraction_aka_numerator);
    mpz_add(the_up_of_the_fraction_aka_numerator, three_x_square, a );

    mpz_t numerator_mul_invert;
    mpz_init(numerator_mul_invert);
    mpz_mul(numerator_mul_invert, the_up_of_the_fraction_aka_numerator, invert);

    mpz_mod(s, numerator_mul_invert,p);

    mpz_t x2; //we want x2 = ((s *s ) -2 * x) % p;
    mpz_init(x2);

    mpz_t s_square;
    mpz_init(s_square);
    mpz_mul(s_square, s, s);

    mpz_t two_x;
    mpz_init(two_x);
    mpz_mul_si(two_x, x, 2);

    mpz_t s_square_minus_two_x;
    mpz_init(s_square_minus_two_x);
    mpz_sub(s_square_minus_two_x, s_square, two_x);

    mpz_mod(x2, s_square_minus_two_x, p);

    mpz_t y2; //we want y2 = (-y + s * (x - x2)) % p;
    mpz_init(y2);

    mpz_t x_minus_x2;
    mpz_init(x_minus_x2);
    mpz_sub(x_minus_x2, x, x2);

    mpz_t s_mul_x_minus_x2;
    mpz_init(s_mul_x_minus_x2);
    mpz_mul(s_mul_x_minus_x2, s, x_minus_x2);

    mpz_t y_plus_s_mul_x_minux_x2;
    mpz_init(y_plus_s_mul_x_minux_x2);
    mpz_t minus_y;
    mpz_init(minus_y);
    mpz_neg(minus_y, y);
    mpz_add(y_plus_s_mul_x_minux_x2, minus_y, s_mul_x_minus_x2);

    mpz_mod(y2, y_plus_s_mul_x_minux_x2, p);

    mpz_set(result.x, x2);

    mpz_set(result.y, y2);

    return result;


}

Coordinates add_point(mpz_t x,mpz_t y,mpz_t x2,mpz_t y2,mpz_t a,mpz_t p){
    Coordinates result;

    mpz_init(result.x);
    mpz_init(result.y);

    if (mpz_cmp_ui(x , 0) == 0 && mpz_cmp_ui(y , 0) == 0){
        mpz_set(result.x, x2);
        mpz_set(result.y, y2);

        return result;

    }

    if (mpz_cmp_ui(x2 , 0) == 0 && mpz_cmp_ui(y2 , 0) == 0){
        mpz_set(result.x, x);
        mpz_set(result.y, y);

        return result;

    }

    if (mpz_cmp(x , x2) == 0 && mpz_cmp(y, y2) == 0) {
        return double_point(x, y, a, p);
    }

    mpz_t invert;
    mpz_t x_minus_x2, x3, y3, s, t;
    
    mpz_init(invert);
    mpz_init(x_minus_x2);
    mpz_init(x3);
    mpz_init(y3);
    mpz_init(s);
    mpz_init(t);

    mpz_sub(x_minus_x2, x, x2);
    mpz_invert(invert, x_minus_x2, p);

    // s = (y1 - y2) * invert % p
    mpz_sub(s, y, y2);
    mpz_mul(s, s, invert);
    mpz_mod(s, s, p);

    // t = (y2 * x1 - y1 * x2) * invert % p
    mpz_mul(t, y2, x);
    mpz_submul(t, y, x2);
    mpz_mul(t, t, invert);
    mpz_mod(t, t, p);

    // x3 = (s^2 - x1 - x2) % p
    mpz_mul(x3, s, s);
    mpz_sub(x3, x3, x);
    mpz_sub(x3, x3, x2);
    mpz_mod(x3, x3, p);

    // y3 = (-s * (s^2 - x1 - x2) - t) % p
    mpz_mul(y3, s, s);
    mpz_sub(y3, y3, x);
    mpz_sub(y3, y3, x2);
    mpz_neg(y3, y3);
    mpz_mul(y3, y3, s);
    mpz_sub(y3, y3, t);
    mpz_mod(y3, y3, p);

    if (mpz_cmp_ui(x3, 0) < 0) {
        mpz_add(x3, x3, p);
    }
    if (mpz_cmp_ui(y3, 0) < 0) {
        mpz_add(y3, y3, p);
    }

    mpz_init_set(result.x, x3);
    mpz_init_set(result.y, y3);

    mpz_clear(invert);
    mpz_clear(x_minus_x2);
    mpz_clear(x3);
    mpz_clear(y3);
    mpz_clear(s);
    mpz_clear(t);

    return result;
}

Coordinates produit_scalaire( mpz_t x, mpz_t y, mpz_t scalar, mpz_t a, mpz_t p) {

    Coordinates temp;
    mpz_init(temp.x);
    mpz_init(temp.y);

    mpz_set_str(temp.x, "0",10);
    mpz_set_str(temp.y, "0",10);

    size_t binary_size = mpz_sizeinbase(scalar, 2);

    char str_scalar_binary[binary_size];
    mpz_get_str(str_scalar_binary, 2, scalar);

    gmp_printf("Scalar int: %Zd\n", scalar);
    
    printf("Scalar binary: %s\n", str_scalar_binary);

    for (size_t i = 0; i < strlen(str_scalar_binary); ++i) {


        temp = double_point(temp.x , temp.y, a, p); 


        if (str_scalar_binary[i] == '1') {
            temp = add_point(temp.x , temp.y , x, y, a, p); 

        }
    }

    return temp;

}


int main() {

    Coordinates result;

    

    mpz_t a;
    mpz_init(a);
    mpz_set_str(a, "0", 10);

    mpz_t scalar;
    mpz_init(scalar);
    mpz_set_str(scalar,"19", 10);

    printf("Parameters of the curve: \n");

    mpz_t x;
    mpz_init(x);
    mpz_set_str(x,"55066263022277343669578718895168534326250603453777594175500187360389116729240",10 );
    gmp_printf("x : %Zd\n", x);

    

    mpz_t y;
    mpz_init(y);
    mpz_set_str(y,"32670510020758816978083085130507043184471273380659243275938904335757337482424",10 );
    gmp_printf("y : %Zd\n", y);

    mpz_t p;
    mpz_init(p);
    mpz_set_str(p,"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
    gmp_printf("p : %Zd\n", p);

    

    result = produit_scalaire(x, y, scalar, a, p);

    print_coordinates(&result);

    return 0;
}
