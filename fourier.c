#include <math.h>

#include "fourier.h"

void normalize(complex s[], int n) {
    for (int k = 0; k < n; k++) {
        s[k].a /= n;
        s[k].b /= n;
    }
}

void nft(complex s[], complex t[], int n, int sign) {
    for (int k = 0; k < n; k++) {
        t[k].a = 0;
        t[k].b = 0;

        for (int j = 0; j < n; j++) {
            double x = sign * 2 * PI * k * j / n;

            double cosx = cos(x);
            double sinx = sin(x);

            t[k].a += s[j].a * cosx - s[j].b * sinx;
            t[k].b += s[j].a * sinx + s[j].b * cosx;
        }
    }
}

void nft_forward(complex s[], complex t[], int n) {
    nft(s, t, n, -1);
}

void nft_inverse(complex t[], complex s[], int n) {
    nft(t, s, n, 1);
    normalize(s, n);
}


//Corrigir (Adicionei operação de sobre escrição para evitar perda de dados devido a uso dos mesmo vetores como entrada e saida)
void fft(complex s[], complex t[], int n, int sign) {
    // fim da recursividade
    if (n == 1) {
        t[0] = s[0];
        return;
    }

    // declarando
    complex sp[n / 2];
    complex si[n / 2];
    complex tp[n / 2];
    complex ti[n / 2];

    // construindo sp e si como foi pedido
    for (int i = 0; i < n / 2; i++) {
        sp[i] = s[2 * i];
        si[i] = s[2 * i + 1];
    }

    //recursividade 
    fft(sp, tp, n / 2, sign);
    fft(si, ti, n / 2, sign);

    for (int k = 0; k < n / 2; k++) {
        
        // mesma formula da nft
        double x = sign * -2.0 * PI * k / n;

        // e^ix = cos(x) +isen(x)
        complex w = {cos(x), sin(x)};
        complex temp;

        // temp usado para armazenar o valor de e^(-2 * pi * ki /n )
        temp.a = w.a * ti[k].a - w.b * ti[k].b;
        temp.b = w.a * ti[k].b + w.b * ti[k].a;

        //t(k ) = tp(k) + temp (que tem o exponencial de e)
        t[k].a = tp[k].a + temp.a;
        t[k].b = tp[k].b + temp.b;
        
        //t(k + n/2) = tp(k) - temp (que tem o exponencial de e)
        t[k + n / 2].a = tp[k].a - temp.a;
        t[k + n / 2].b = tp[k].b - temp.b;
        
    }

}

void fft_forward(complex s[], complex t[], int n) {
    fft(s, t, n, -1);
}

void fft_inverse(complex t[], complex s[], int n) {
    fft(t, s, n, 1);
    normalize(s, n);
}

void fft_forward_2d(complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
    // Vetor auxiliar 
    complex temp[MAX_SIZE];

    // FFT nas linhas primeiro
    for (int y = 0; y < height; y++) {
        fft_forward(matrix[y], temp, width);
        for (int x = 0; x < width; x++) {
            matrix[y][x] = temp[x]; // Sobrescreve a entrada
        }
    }

    // FFT nas colunas depois
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            temp[y] = matrix[y][x]; 
        }
        fft_forward(temp, temp, height); 
        for (int y = 0; y < height; y++) {
            matrix[y][x] = temp[y]; // Sobrescreve a entrada
        }
    }

}

//Verificar se está funcionando (inverti a ordem dos loops de linhas e colunas)
//Vetor auxiliar so vai até maior dimensão da matrix, como ele é sobre escrito a cada interação dos loops que o usam
void fft_inverse_2d(complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
    //Exatamente o mesmo conceito da fft_forward_2d porem muda a funçao chamada pela fft_inverse
    
    // Vetor auxiliar
    int size = 0;
    if(width > height){
        size = width;
    }
    else{
        size = height;
    }

    complex temp[size];

    // FFT nas colunas 
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            temp[y] = matrix[y][x]; 
        }
        fft_inverse(temp, temp, height); 
        for (int y = 0; y < height; y++) {
            matrix[y][x] = temp[y]; // Sobrescreve a entrada
        }
    }

    // FFT nas linhas 
    for (int y = 0; y < height; y++) {
        fft_inverse(matrix[y], temp, width); 
        for (int x = 0; x < width; x++) {
            matrix[y][x] = temp[x]; // Sobrescreve a entrada
        }
    }

}

void filter(complex input[MAX_SIZE][MAX_SIZE], complex output[MAX_SIZE][MAX_SIZE], int width, int height, int flip) {
    int center_x = width / 2;
    int center_y = height / 2;

    double variance = -2 * SIGMA * SIGMA;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int dx = center_x - (x + center_x) % width;
            int dy = center_y - (y + center_y) % height;

            double d = dx * dx + dy * dy;
            double g = exp(d / variance);

            if (flip) {
                g = 1 - g;
            }

            output[y][x].a = g * input[y][x].a;
            output[y][x].b = g * input[y][x].b;
        }
    }
}

void filter_lp(complex input[MAX_SIZE][MAX_SIZE], complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
    filter(input, output, width, height, 0);
}

void filter_hp(complex input[MAX_SIZE][MAX_SIZE], complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
    filter(input, output, width, height, 1);
}
