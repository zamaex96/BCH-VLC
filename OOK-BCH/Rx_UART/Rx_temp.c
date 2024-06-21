
#define F_CPU 16000000
#define CHECKBIT(x,y) (x & (y))

#include <stdio.h>
#include <avr/io.h>
//#include <mega128.h>
#include <util/delay.h>
#include <stdlib.h>
#include <math.h>		// For use of power function

//** Note: Must be compiled with "char is unsigned" **//

//** Variable Definitions **


//Global Vars
unsigned long ReceivedMessage;	// Encoded Message as a long
char EncMsgArray[32];			// Encoded Message in Array format
char n=31, k=16, t=3;			// BCH Code definitions
unsigned long GeneratorPoly = 36783;	// Generator Polynomial: 107657 in octal

// Struct for use in Berlekamp Decoding, allows for arrays of Poly32's
// Note that all GF32 arrays contain the POWERS of alpha as their data,
// not the actual GF32 element - hence multiplication is actually addition
struct Poly32
{
	char p[32];
};


// Define GF32 - ZERO as negative one since a zero would
// zero out various multiplications
#define ZERO -1

// Efficient lookup tables for GF32
// Alpha exponents - computed as powers of 2
// Where alpha is from the Galois Field(32)
 char lookup[] = {
                1,  // a^0  = 0
                2,  // a^1  = 1
                4,  // a^2  = 2
                8,  // a^3  = 3
                16, // a^4  = 4
                5,  // a^5  = 2 0
                10, // a^6  = 3 1
                20, // a^7  = 4 2
                13, // a^8  = 3 2 0
                26, // a^9  = 4 3 1
                17, // a^10 = 4 0
                7,  // a^11 = 2 1 0
                14, // a^12 = 3 2 1
                28, // a^13 = 4  3 2
                29, // a^14 = 4 3 2 0
                31, // a^15 = 4 3 2 1 0
                27, // a^16 = 4 3 1 0
                19, // a^17 = 4  1 0
                3,  // a^18 = 1 0
                6,  // a^19 = 2 1
                12, // a^20 = 3 2
                24, // a^21 = 4 3
                21, // a^22 = 4 2 0
                15, // a^23 = 3 2 1 0
                30, // a^24 = 4 3 2 1
                25, // a^25 = 4 3 0
                23, // a^26 = 4 2 1 0
                11, // a^27 = 3 1 0
                22, // a^28 = 4 2 1
                9,  // a^29 = 3 0
                18, // a^30 = 4 1
                1   // a^31 = 0
            };

 char reverseLookup[] = {
                ZERO, // a^0
                0,  // a^1
                1,  // a^2
                18, // a^3
                2,  // a^4
                5,  // a^5
                19, // a^6
                11, // a^7
                3,  // a^8
                29, // a^9
                6,  // a^10
                27, // a^11
                20, // a^12
                8,  // a^13
                12, // a^14
                23, // a^15
                4,  // a^16
                10, // a^17
                30, // a^18
                17, // a^19
                7,  // a^20
                22, // a^21
                28, // a^22
                26, // a^23
                21, // a^24
                25, // a^25
                9,  // a^26
                16, // a^27
                13, // a^28
                14, // a^29
                24, // a^30
                15  // a^31
            };

/** Function Prototypes **/
// Note: For GF2 unsigned longs can be used since only two elements
// in the field.  For GF32, byte arrays needed to be used since there
// are 32 elements.


// Initialization function
void Initialize();

// Main Encoding/Decoding System
void System();

// All encoding for BCH ECC is done here
void DecoderBCH();

/** Helper Functions **/

// Finds the degree of a polynomial encoded as a long in GaloisField2
unsigned char GF2FindDegree(unsigned long a);

// Adds 2 polynomials encoded as a long - GF2
unsigned long GF2Add(unsigned long a, unsigned long b);

// Polynomial Multiplication for longs in GF2
unsigned long GF2Multiply(unsigned long a, unsigned long b);

// Polynomial Division for longs in GF2
void GF2Divide(unsigned long a, unsigned long b, unsigned long *qr);

// Retrieves a specified bit from a long
unsigned char getBit(unsigned long r, char i);

// Converts a long into a 32 length byte array
void Bits2Bytes(unsigned long num, char *p);

// Initializes GF32 arrays to ZERO
void GF32Init(char *p);

// Adds two alpha coefficients in GF32
char GF32add2alpha(unsigned char a, unsigned char b);

// Finds the degree of a GF32 polynomial
char GF32FindDegree(char *p);

// Evaluates the result of a GF32 Polynomial for some x
char GF32Evaluate(char a, char *p);

// Adds all alpha coefficients pairwise in 2 Arrays in GF32
void GF32Add(char *a, char *b, struct Poly32 *powers);

// Polynomial Multiplication for longs in GF32
void GF32Multiply(char *a, char *b, struct Poly32 *mul);

// Multiplies a GF32 polynomial with some power of x
void multiplyX(char x_power, char *p, struct Poly32 *ret);

// Multiplies a GF32 polynomial by a constant
void multiplyConstant(char c, char *p, struct Poly32 *ret);

// Corrects the detected errors in an encoded message
void CorrectErrors(char *p);

// Converts an array of bytes into a long
unsigned long Bytes2Bits(char *p);

// Parses the lower 16bits of a long into 2 chars
void deConcate(unsigned long a, unsigned char *ret);

// Prints a long as 4 chars
void TransLong(unsigned long a);

void transmit(unsigned char sentdata);



//**Initialize**
//Sets up vars, timers, and Mega32 registers
void Initialize()
{

	n=31; k=16; t=3;	// Initialize BCH Code parameters to (31,16)
	//ReceivedMessage=0;
	GeneratorPoly = 36783; // Init the Generator polynomial
	//-------------------------------------------------------------

	unsigned int baudrate[15] ={416,207,103,68,51,34,25,16,12,8,7,3,3,1,0};
 
 // UBRR1L = baudrate[2];		//set baudrate - cast low byte
 // UBRR1H = baudrate[2]>>8;	//set baudrate - cast high byte 9600bps  
  //UBRR1H = baudrate[3]>>8;				//11400bps

   
  
  UBRR1L = baudrate[12];
  UBRR1H = baudrate[12]>>8;               // 4:19200bps   //7:57600 //8:76800 // 9:115200

  UBRR0L = baudrate[12];
  UBRR0H = baudrate[12]>>8;               //19200bps  
   
  //UBRR=((f_clk/(BaudRate*16))-1)
 
  
  // UCSR1A = (1 << U2X0);                //double baudrate 
  UCSR1C = (1 << UCSZ1) | (1 << UCSZ0);  //serial 8-bit format. no parity, stop bit 1, data 8
  UCSR1B = (1 << TXEN1) | (1 << RXEN1);    //enable tx rx located in control reg UCSR1B
                                         //When set, these two bits turn on the serial buffers to allow for serial communications
//UCSR1A = (1 << U2X1);                //double baudrate 
   //  UCSR0A = (1 << U2X0);
  UCSR0C = (1 << UCSZ1) | (1 << UCSZ0);
  UCSR0B = (1 << TXEN0) | (1 << RXEN0);
  // #if USE_2X
     UCSR0A |= (1 << U2X0); //double baudrate
	 UCSR1A |= (1 << U2X1); //double baudrate
     //#else
     //UCSR0A &= ~(1 << U2X0)
}

//**Encode/Decode System**
// High level system for the encoding and decoding operations for
// the entire message to be transmitted
void System()
{
	 unsigned char out[2];
	unsigned long decoded_bits;	// Final decoded message of 2 chars
    
		GF32Init(EncMsgArray); // Initialize to ZERO's
	Bits2Bytes(ReceivedMessage, EncMsgArray);	// Convert to a GF32 polynomial
	    
		DecoderBCH();	// and Decoded the resulting encoded message
		decoded_bits = Bytes2Bits(EncMsgArray);
		deConcate(decoded_bits, out);	// Parse out the 2 chars

			transmit(out[0]);
			transmit(out[1]);

}

//** 32 bit BerleKamp Decoder for BCH Codes **
// Decodes an encoded 32 bit message, and corrects for any detected errors
// according to the Berlekamp algorithm. See the writeup for details.
void DecoderBCH()
{
	unsigned char i, j;
	char k=0;
	char delta[7];
	char Syndromes[32];
 	//Lambda is an array of 7 GF32 polynomials structs
	//Note: Lambda/T arrays can be small (i.e length 4) if rewrite multiply/divide, etc.
	struct Poly32 lambda[7];
	struct Poly32 T[7];
	struct Poly32 temp;

	GF32Init(Syndromes);
	// Create Syndrome Polynomial
	for (i=0; i<2*t; i++)
		Syndromes[i+1] = (char)GF32Evaluate(i+1, EncMsgArray);


   /* 1  - Initialization */
   	// add 1 to S(x) and initialize Berlekamp's Algorithm
	Syndromes[0] = 0;

	//Init Lambda[i] polynomials
	for (i=0; i<7; i++)
		for (j=0; j<32; j++)
		{
			lambda[i].p[j]=ZERO;
			T[i].p[j]=ZERO;
		}

	// lambda_0 (x) = 1
	lambda[0].p[0] = 0;
    	// T_0 (x) = 1
    	T[0].p[0] = 0;

	while( k < t )
	{
      /* Berlekamp Algorithm */

	  /* 2 */    // Delta[2k] = coeff. of x^(2k+1) in the product Lambda[2k](x) * [1 + Syn(x)]
			 GF32Multiply(lambda[2*k].p, Syndromes, &temp);
			 delta[2*k]  = (char)temp.p[2*k+1];

	  /* 3 */ 	 // Lambda[2k+2](x) = Lambda[2k](x) + Delta[2k]*(x*T[2k](x))
			 multiplyX(1, T[2*k].p, &temp);
			 multiplyConstant(delta[2*k], temp.p, &temp);
			 GF32Add(lambda[2*k].p, temp.p, &lambda[2*k+2]);

	  /* 4 */
			 if (delta[2*k] == ZERO || (char)GF32FindDegree(lambda[2*k].p) > k)
				multiplyX(2, T[2*k].p, &T[2*k+2]);
			 else
			 {
				multiplyX(1, lambda[2*k].p, &temp);
				multiplyConstant((char)(31-delta[2*k]), temp.p, &T[2*k+2]);
			 }


	  /* 5 */     k++; // Increment for next iteration
        }

	  CorrectErrors(lambda[2*k].p);	// Correct the errors as determined by the locater
}							// Lambda polynomial

//** Find Degree of a Polynomial in GF2 **//
// Simply finds the first index of a 32bit number that is not zero
// And that is the degree+1 of the polynomial
unsigned char GF2FindDegree(unsigned long num)
{
	unsigned char i=0, deg=0;

	num = num << 1;	// Shift left since top bit is ignored in algorithm
	for(i=0; i<30; i++)
	{
		if (num & 0x80000000)	// Mask the current top bit, to see if it's a one
			{deg= (30-i); break;}	// if so, that's the degree
		num = num << 1;		// otherwise, keep shifting
	}
	return deg;
}

//** Polynomial Addition in GF2 **//
// Simply executes Modulo 2 addition
unsigned long GF2Add(unsigned long a, unsigned long b)
{
        return (a^b);	// simply xor the bits (GF2 addition for polynomials)
}

//** Polynomial Multiplication in GF2 **//
// Executes Multiplication in GF2 for polynomials
unsigned long GF2Multiply(unsigned long a, unsigned long b)
{
    unsigned long mul = 0;
    unsigned long add;

	char i;

	add = b;

 	for(i=0; i <= GF2FindDegree(a); i++) // loop while not to the end of the poly
	{
		if(getBit(a, i) == 1)		// If coeff. is a one, then add multiplicand
			mul ^= add;
		add = add<<1;			// and shift the multiplicand up one
  	}

    return mul;
}

//** Polynomial Long Division in GF2 **//
// Executes Long Division in GF2 for polynomials
// The remainder (qr[1]) is equal to the final dividend (qr[0])
// Degree of qr[1] should be smaller than degree of divisor (the break condition in the loop)
void GF2Divide(unsigned long a, unsigned long b, unsigned long *qr)
{
	unsigned long dividend;
	unsigned long divisor;

	unsigned long q;
 	int deg = 0;

	dividend = a;
	divisor = b;
	qr[0] = 0;

	while(1)	// Keep doing this until break is activated
	{
		// Subtract degrees to find what the degree of each term in the quotient
		deg = (int)(GF2FindDegree(dividend) - (int)GF2FindDegree(divisor));

		if (deg < 0) 	// If negative, then you are done
		{
			qr[1] = dividend;		// return the dividend as the remainder
			return ;
		}

		if (deg > 0)	// otherwise find the appropriate degree for the term
			q = (unsigned long)pow((float)2,(float)deg)+1;
		else
			q = 1;
		qr[0] = GF2Add(qr[0], q);	// and add the term to the quotient
		// finally, reduce (i.e add mod 2) the divided by (term*divisor)
		dividend = GF2Add(dividend, (GF2Multiply(q, divisor)));
	}
	qr[1] = dividend;		// Return the remainder
}

//** Get a bit from a Long **//
// Returns the bit i of the long r
unsigned char getBit(unsigned long r, char i)
{
        unsigned char ret;

	  // Shifts and Masks to get the appropriate bit
        ret = ((r<<(32-i-1))>>31)& 0x00000001;
        return ret;
}

//** Long to Array convertor **//
// Takes a polynomial in GF2 (long) and coverts it into a polynomial in
// GF32 (a byte array)
void Bits2Bytes(unsigned long num, char *p)
{
	unsigned char i=0, temp=0;

 	for(i=0; i<32; i++)
	{
		temp = num % 2;
		if (temp == 0)
			p[i] = ZERO;	// -1 is ZERO, i.e. coeff = 0
		else
			p[i] = 0;		// alpha**0, i.e. coeff = 1
		num = num >> 1;		// shift for next iteration
	}
}

//** GF32 Initialize **//
// Simply initializes a GF32 array to all ZERO's
void GF32Init(char *p)
{
	 unsigned char i=0;

	 for (i=0; i<32; i++)
		 p[i]=ZERO;
}




//** Add Two Alpha Coeff. in GF32 **//
// Uses the precomputed lookup tables to add powers of alpha mod 32
char GF32add2alpha(unsigned char a, unsigned char b)
{
        if ((a == ZERO) && (b == ZERO))	// ZERO+ZERO=ZERO
            return ZERO;
        if (a == ZERO)				// ZERO is additive identity
            return b;
        else if (b == ZERO)
            return a;
        else					// Simply XOR and use lookup
            return reverseLookup[lookup[a]^lookup[b]];
}

//** Find Degree of GF32 Polynomial **//
// Returns the index of the first non-ZERO element
char GF32FindDegree(char *p)
{
       unsigned char i = 32;
        while(--i > 0)
            if (p[i] != ZERO) return i;
        return 0;
}

//** GF32 Polynomial Evaluation **//
// Evaluates the result of a polynomial defined over GF32 evaluated
// with an element from GF32
char GF32Evaluate(char a, char *p)
{
	char ret = ZERO; unsigned char i=0;
	char pow=0;

      for(i=0; i <= GF32FindDegree(p); i++) // evaluate over the length of the polynomial
	{
          if (p[i] != ZERO)
	    {
                pow = (char)((p[i]+a*i) % 31); // index is the degree, multiply exponents
                if (pow < 0)
                	pow = (char)(31+pow); // Evaluate mod 32
                ret = GF32add2alpha(ret, pow);	// exponent multiplication = add
            }
        }
        return ret;
 }

//** GF32 Add Two Polynomials **//
// Adds two GF32 polys using the lookup tables for each pairwise coeff.
void GF32Add(char *a, char *b, struct Poly32 *powers)
{
	unsigned char i=0;

        for (i=0; i < 32; i++)
            powers->p[i] = GF32add2alpha(b[i], a[i]);
}

//** GF32 Polynomial Multiplication **//
// Multiplies two GF32 polynomials and returns the result by reference.
void GF32Multiply(char *a, char *b, struct Poly32 *mul)
{
	struct Poly32 add;
	unsigned char i,j;

	for(i=0; i<32; i++)	// Initialize the arrays
	{
		mul->p[i]=ZERO;
		add.p[i]=ZERO;
	}

      for(i=0; i <= GF32FindDegree(a); i++)
	{
          if(a[i] != ZERO )	// multiply only non-zero terms
	    {
            for(j=0; j <= GF32FindDegree(b); j++) // add then shift
		{
                    if(b[j] != ZERO)
                        add.p[j+i] = (char)((a[i]+b[j]) % 31);
            }
            GF32Add(mul->p, add.p, mul);
            GF32Init(add.p);
          }
      }
}

//** Multiply GF32 Polynomial by x^power **//
// Simply executes a cyclic shift (mod 32) by the power of x
void multiplyX(char x_power, char *p, struct Poly32 *ret)
{
	unsigned char i;

        for(i=0; i<32; i++)
            ret->p[(i+x_power) % 32] = p[i];	// cyclic shift mod 32
}

//** Multiply GF32 Polynomial by Constant **//
// Simply multiplies each coeff. by an element from GF32, mod 32
void multiplyConstant(char c, char *p, struct Poly32 *ret)
{
	 unsigned char i;

      // if multiplying by zero, return zero
      if(c == ZERO)
	{
		GF32Init(ret->p);
		return ;
	}

      for (i = 0; i < 32; i++)
	{
            if(p[i] != ZERO )
                ret->p[i] = (char)((p[i]+c) % 31); // add the constant exponent, mod 32
            else
                ret->p[i] = ZERO;
      }
}

//** Correct the Errors in the Encoded Message **//
// Corrects the encoded message according to the errors detected by the locater
// polynomial, Lambda which is passed as p - i.e. the roots of Lambda
// are the locations of the errors
void CorrectErrors(char *p)
{
	unsigned char i;

	// evaluate roots of lambda[2*k] and flip the received code word bits accordingly
	for(i=0; i<32; i++)
	{
     		if (GF32Evaluate(i, p) == ZERO )	// Find the roots of Lambda
	  	{
       		
         		if (EncMsgArray[31-i] == ZERO)	// Simply flip the bits
         			EncMsgArray[31-i] = (char)0;
         		else
         			EncMsgArray[31-i] = ZERO;
        }
   	}
}

//** Convert GF32 polynomial (array) to a GF2 polynomial (long) **//
// Simply Or-in the appropriate bits, and shift up
unsigned long Bytes2Bits(char *p)
{
	char i;
	unsigned long ret=0;

	for(i=0; i<31; i++)
	{
		ret = ret | (p[31-i] == 0); // if 0, or in a 1, if ZERO, or in a 0
		ret = ret << 1;	// and then shift it up
	}
	ret = ret | (p[0] == 0);	// shift up the last one - since only 31 elements

	return ret;
}

//** De-Concatenate **//
// Returns the bottom 16 bits of a long as 2 chars
void deConcate(unsigned long a, unsigned char *ret)
{
	a = a << 1;
	ret[0] = (unsigned char)((a & 0xFF000000) >> 24);
	ret[1] = (unsigned char)((a & 0x00FF0000) >> 16);
}





//-----------------------------------------------------------------------------------------------------------

/*
void adc_init()
{
	ADCSRA=0X00;
	ADMUX = (1<<REFS0);
	ADCSRA= (1<<ADEN | 1<<ADPS0 | 1<<ADSC | 1<<ADFR );
	//ADCSRA= 0xC7;
}
*/

void transmit (unsigned char sentdata)
{
	// Wait for empty transmit buffer 
	while ( !(UCSR1A & (1<<UDRE1)) );

	// Put data into buffer, sends the data 
	UDR1 = sentdata;
	
	
}

unsigned char Rxdata (void)
{

	// Wait for empty transmit buffer 
	while ( !(UCSR0A & (1<<RXC0)) );

	// Receive data from buffer, return the data 
	return UDR0;	

}

/*
//subroutine to detect the first 8 ones
unsigned int detect_all_one (void)
{

	unsigned int temp; //to store PINA info
	unsigned int byte_placeholder; //the bucket to hold 8 bits of PINA info
	unsigned int bit_to_shiftin; //shift in 1 or 0 depending on PINA
	unsigned int bit_to_check; //check the last bit of PINA
	bit_to_check = 0x02;
	byte_placeholder = 0x00;

	while(byte_placeholder != 0xff)
  	{

		temp = PINA;
		if (!(CHECKBIT(temp, bit_to_check)))	// Note that PINx is a inverted signal
		{
			bit_to_shiftin = 0x00;
			byte_placeholder = byte_placeholder << 1; //shift the placeholder to the left by one bit
			byte_placeholder = byte_placeholder | bit_to_shiftin;
			//PORTF = 0x01;
			//_delay_ms(300);
			//PORTF = 0x02;
			//_delay_ms(300);
 		}
		else
		{
			bit_to_shiftin = 0x01;
			byte_placeholder = byte_placeholder << 1; //shift the placeholder to the left by one bit
			byte_placeholder = byte_placeholder | bit_to_shiftin;
			
		}
		_delay_us(500);
		
  	}

	return 1;
}
*/

unsigned int detect_all_one (void)
{
	unsigned char byte_data;
	byte_data = 0x00;

	while(byte_data !=0xff)
	{
		byte_data = Rxdata();		
	}

	return 1;
	
}





//subroutine for start of packet detection
unsigned int detect_startofpacket (void)
{
	unsigned int second_sequence; //the status flag to ensure 10000000 are detected
	//second_sequence = Rxdata();
	second_sequence = 0x00;
	while(second_sequence != 0x80)
	{
		detect_all_one();
		second_sequence = Rxdata();

		
	}
	return 1;
}

//** Print a long as 4 chars **//
// For use with Hyperterm, since it cannot display longs
void TransLong(unsigned long a)
{
	unsigned char b,c,d,e;

	b = (unsigned char)((a & 0xFF000000) >> 24);
	c = (unsigned char)((a & 0x00FF0000) >> 16);
	d = (unsigned char)((a & 0x0000FF00) >> 8);	
	e = (unsigned char)((a & 0x000000FF));
	transmit(0xFF);
	//transmit_ASK(0xFF);
	transmit(0x80);
	//transmit_ASK(0x80);
	transmit(b);
	//transmit_ASK(b);
	transmit(c);
	//transmit_ASK(c);
	transmit(d);
	//transmit_ASK(d);
	transmit(e);
//	transmit_ASK(e);
	//transmit(0x00);
	//transmit_ASK(0x00);
}

int main(void)
{
	//unsigned int byte_placeholder;
	//unsigned int key;
	unsigned int ctr,i;
	unsigned char data[5];
	unsigned long Det_long[2],data_long[5];
          
	DDRA = 0xFF;
	PORTA = 0x01;

	Initialize();

	//PORTB = 0xF0;
	
	while(1)
	{
		
		for(i=0; i<2;i++)
		{
		    detect_startofpacket();
		
		    //byte_placeholder = shift_byte();
		    //transmit(0xFF); 	// Start byte
		    //transmit(0x80);// Header byte
		    for(ctr=0; ctr<4;ctr++)
		    {
			  data[ctr] = Rxdata();
			  //transmit(data[ctr]);
			  data_long[ctr]=data[ctr];
		    }
		   //transmit(0x00);
           //Det_long= ((data[0] & 0x000000FF) << 24)|((data[1] & 0x000000FF) << 16)|((data[2] & 0x000000FF) << 8)| ((data[3] & 0x000000FF));
            Det_long[i]= (data_long[0]<<24) | (data_long[1] << 16) | (data_long[2] << 8) | (data_long[3]);
		   //TransLong(Det_long[i]);
		   
		   ReceivedMessage=Det_long[i];
		   TransLong(ReceivedMessage);
		  // System();
		   
	    }
			

		

	} 
	return 0;
}





