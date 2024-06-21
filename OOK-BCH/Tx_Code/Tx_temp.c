
#define F_CPU 16000000
#define CHECKBIT(x,y) (x & (y))
#include <stdio.h>
#include <avr/io.h>
#include <avr/interrupt.h>
#include <stdlib.h>
#include <math.h>
#include <util/delay.h>

volatile long cnt;
char Message[] = "AI\r\n";
//Global Vars
unsigned long EncodedMessage=0;	// Encoded Message as a long
//unsigned long long_placeholder=0xFFFFFFFF; 
char EncMsgArray[32];			// Encoded Message in Array format
char Enc_message[32];
char n=31, k=16, t=3;			// BCH Code definitions
unsigned long GeneratorPoly = 36783;	// Generator Polynomial: 107657 in octal
char seed=0;				// Random seed for Error module
unsigned long qr[2] = {0,0};
//Counters to keep track of errors for BER calculations
int trans_error_count=0,rec_error_count=0;
int enc_message_bits;
// Struct for use in Berlekamp Decoding, allows for arrays of Poly32's
// Note that all GF32 arrays contain the POWERS of alpha as their data,
// not the actual GF32 element - hence multiplication is actually addition
struct Poly32
{
	char p[32];
};
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
// Note: For GF2 unisgned longs can be used since only two elements
// in the field.  For GF32, byte arrays needed to be used since there
// are 32 elements.
// Main Encoding/Decoding System
void System();
// All encoding for BCH ECC is done here - accepts 2 char message
void EncoderBCH(unsigned char a, unsigned char b);
// All decoding for BCH ECC is done here
void DecoderBCH();
/** Helper Functions **/

// Concatenates 2 chars into a long
unsigned long Concate(unsigned char num1, unsigned char num2);

// Finds the degree of a polynomial encoded as a long in GaloisField2
unsigned char GF2FindDegree(unsigned long a);

// Adds 2 polynomials encoded as a long - GF2
unsigned long GF2Add(unsigned long a, unsigned long b);

// Polynomial Multiplication for longs in GF2
unsigned long GF2Multiply(unsigned long a, unsigned long b);

// Polynomial Division for longs in GF2
void GF2Divide(unsigned long a, unsigned long b, unsigned long *qr);

// Polynomial Modulo for longs in GF2
unsigned long GF2Mod(unsigned long a, unsigned long b);

// Retrieves a specified bit from a long
unsigned char getBit(unsigned long r, char i);

// Converts a long into a 32 length byte array
void Bits2Bytes(unsigned long num, char *p);

// Initializes GF32 arrays to ZERO
void GF32Init(char *p);

// Prints an Array byte by byte
void TransArray(char *p);

// Prints an array in GF32 polynomial format
//void GF32PrintArray(char *p);

// Adds two alpha coefficients in GF32
char GF32add2alpha(char a, char b);

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

void delay_func (unsigned int over);
// Prints a long as 4 chars
void TransLong(unsigned long a);
//transmit 32 bit long via USART
void transmit0(unsigned char data);
void transmit1(unsigned char data);
void transmit_ASK (unsigned char storage);
void transmit_32b (unsigned long storage);
 //Initialization function
void Initialize();

//**Initialize**
//Sets up vars, timers, and Mega32 registers
void Initialize()
{

	n=31; k=16; t=3;	// Initialize BCH Code parameters to (31,16)
	EncodedMessage=0;
	GeneratorPoly = 36783; // Init the Generator polynomial
			
	

}

void System()
{ 	
	unsigned char a, b;
	unsigned int CurrIndex = 0;	// Keeps track of where we are in the Message


	a = Message[CurrIndex];	// Initialize the first 2 chars to be encoded
	b = Message[CurrIndex+1];

	while ((a != '\r') && (b != '\r'))	// Loop throught the message until terminator chars are reached
	{
	
		EncoderBCH(a,b);	// Run the encoder on the first 2 chars of Message
	 	CurrIndex += 2;		// Increment to grab the next 2 chars to be Encoded/Decoded
		a = Message[CurrIndex];
		b = Message[CurrIndex+1];
		enc_message_bits += 31;		// Every iteration, 31 bits are encoded
	}
	
}
void EncoderBCH(unsigned char a, unsigned char b)
{
	// Systematic encoding as follows:
	// (m(x) * x^(n-k)) mod g(x) + (m(x) * x^(n-k))
	
	// Shift by x^(n-k) = x*15
	EncodedMessage = Concate(a, b);	// Concatenates the two chars to a 16bit message
	EncodedMessage = EncodedMessage << (n-k);	// Multiply by x^(n-k)
	
	EncodedMessage = GF2Add(GF2Mod(EncodedMessage, GeneratorPoly), EncodedMessage);

//	EncodedMessage ^= 0x20010008;	// Errors can be manually added to message for debugging
	
	GF32Init(EncMsgArray); // Initialize to ZERO's
	Bits2Bytes(EncodedMessage, EncMsgArray);	// Convert to a GF32 polynomial

	TransArray(EncMsgArray);

}
//** Concatenate **//
// Concatenates two 8bit numbers into a 16bit message
// Since it's a (31,16) code
unsigned long Concate(unsigned char num1, unsigned char num2)
{
	unsigned long temp=0;
	
	temp = temp | num1;	// Or in the num1 and shift it up
	temp = temp << 8;
 	temp = temp | num2;	// Or in num2
	
	return temp;
}

//** Find Degree of a Polynomial in GF2 **//
// Simply finds the first index of a 32bit number that is not zero
// And that is the degree+1 of the polynomial
unsigned char GF2FindDegree(unsigned long num)
{
	char i=0; unsigned char deg=0;
	
	num = num << 1;	// Shift left since top bit is ignored in algorithm
	for(i=0; i<30; i++)
	{
		if (num & 0x80000000)	// Mask the current top bit, to see if it's a one
			{deg=(30-i);break;}	// if so, that's the degree
        else
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

//** Polynomial Modulo in GF2 **//
// Simply executes GF2Divide to find remainder of two polynomials
unsigned long GF2Mod(unsigned long a, unsigned long b)
{
	//unsigned long qr[2] = {0,0};
	GF2Divide(a, b, &qr[0]);
	return qr[1];
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
	int i=0, temp=0;
	
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
	 int i=0;
	 
	 for (i=0; i<32; i++)
		 p[i]=ZERO;
}

//** GF32 Initialize **//
// Prints a 32-element array

void TransArray(char *p)
{
	char i=0;int m=0;//	unsigned int bit_to_shiftin;
	
	for(i=0; i<32; i++)
	{
		//if (i%8 == 0) printf("  ");	// Space every 8 elements for clarity
		if (p[31-i] == ZERO)
		   {	// i.e. if ZERO, send 0
			Enc_message[m]=0;m++;
		//	PORTA = 0x00;
		//	delay_func (5);
			//	bit_to_shiftin = 0x00000000;
				//long_placeholder = long_placeholder << 1;
				//long_placeholder = long_placeholder | bit_to_shiftin;
				//long_placeholder = long_placeholder & 0xFFFFFFFF;
			}
		else
			{
			Enc_message[m]=1;m++;
		//	PORTA = 0x03;
		//	delay_func (5);
				//bit_to_shiftin = 0x00000001;
				//long_placeholder = long_placeholder << 1;
			//	long_placeholder = long_placeholder | bit_to_shiftin;
			//	long_placeholder = long_placeholder & 0xFFFFFFFF;
			}	// otherwise it's a one
	}
TransLong(EncodedMessage);
//	transmit_ASK(0xFF);
//	transmit_ASK(0x80);	
//transmit_32b(EncodedMessage);
//transmit_ASK(0x00);
}


void delay_func (unsigned int overfl)
{
         int cnt1=0;
        cnt1 = 0, cnt=0;
		while(cnt1==overfl)
		{
			cnt1 = cnt;
		}
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
	transmit1(0xFF);
	//transmit_ASK(0xFF);
	transmit1(0x80);
	//transmit_ASK(0x80);
	transmit1(b);
//	transmit_ASK(b);
	transmit1(c);
//	transmit_ASK(c);
	transmit1(d);
//	transmit_ASK(d);
	transmit1(e);
//	transmit_ASK(e);
//	transmit1(0x00);
	//transmit_ASK(0x00);
}







//-----------------------------------serial----------------------------------------
void Timer1_init(void) 
{
  
  	unsigned int baudrate[15] ={416,207,103,68,51,34,25,16,12,8,7,3,3,1,0};
	//unsigned int cycles;
	unsigned long cycles;
	 // unsigned int tcnt1;
  	// Initializing Timer
	TCCR1A = 0;                 // clear control register A 
	TCCR1B |= (1 << WGM13); // set mode 8: phase and frequency correct pwm, stop the timer
 	//do {
	unsigned char oldSREG = SREG;
	cli();
	//cycles = 0xffd0;	//	65536 - 48 = 65488 => 0xffd0
	//TCNT1 = cycles;
	// Top=ICR1=16000000/(2*prescaler*desired interrupt frequency), compare match register = [ 16,000,000Hz/ (prescaler * desired interrupt frequency) ] - 1
	//  ! Remember ! that when you use timers 0 and 2 this number must be less than 256, and less than 65536 for timer1
    // ADJUSTABLE VARIABLES
// Strobe frequency
//uint16_t timer1Prescaler = 64; /* 1, 8, 64, 256, 1024 */
//uint8_t strobePeriod = 50, /* milliseconds */
//  strobeDutyCycle = 20; /* percent */
// Set PWM frequency/top value
//  ICR1 = (F_CPU*strobePeriod / (timer1Prescaler*1000) ) - 1; /* equals 12499 */
//OCR1A = ICR1 / (100 / strobeDutyCycle); /* equals 2499 */
	cycles = 5;	 // 2/(16000000/8)*200 = 0.2ms
	ICR1 = cycles;
	//sei();
	//tcnt1 = TCNT1;
	SREG = oldSREG;
	//}while(tcnt1==0);
    //cycles = 125;	 // 2/(16000000/8)*200 = 0.2ms
	//ICR1 = cycles;
	
	TCCR1B &= ~((1 << CS10) | (1 << CS11) | (1 << CS12));
	TCCR1B |= (1 <<CS11 );	// prescale: clk / 8
	/* Timer clock = I/O clock / 64 */
    //TCCR1B = (1<<CS11)|(1<<CS10);

	UBRR0L = baudrate[12];
    UBRR0H = baudrate[12]>>8;  
	UCSR0C = (1 << UCSZ1) | (1 << UCSZ0);  //serial 8-bit format. no parity, stop bit 1, data 8
    UCSR0B = (1 << TXEN0) | (1 << RXEN0);

	  UBRR1L = baudrate[12]; 
  	  UBRR1H = baudrate[12]>>8;  
	  UCSR1C = (1 << UCSZ1) | (1 << UCSZ0);  // no parity ,  stop bit 1 , data 8
 	  UCSR1B = (1 << TXEN1) | (1 << RXEN1);
	 // #if USE_2X
     UCSR0A |= (1 << U2X0); //double baudrate
	UCSR1A |= (1 << U2X1); //double baudrate
     //#else
     //UCSR0A &= ~(1 << U2X0);

	
}

void transmit0(unsigned char data)
{
while(!(UCSR0A & (1<<UDRE0)));
UDR0=data;

}

void transmit1(unsigned char data)
{
while(!(UCSR1A & (1<<UDRE1)));
UDR1=data;

}

unsigned int Rxdata (void)
{

	// Wait for empty  buffer 
	while ( !(UCSR0A & (1<<RXC0)) );

	// Receive data from buffer, return the data 
	return UDR0;	

}

ISR(TIMER1_OVF_vect)
{
	cnt++;
}

void transmit_ASK (unsigned char storage)	// subroutine to transmit SRAM data to PORTA
{
	int counter = 0;
	unsigned char bit_to_check = 0x80;  //check the right most bit
	//TIMSK = (1 << TOIE1);		// Timer1 Interrupt Mask Enable
	

	while (counter < 8) 
	{
		if (!(CHECKBIT(storage, bit_to_check)))
		{
			PORTA = 0xFE;
      	}
		else
		{
			PORTA = 0xFF;
		}

		delay_func(0);
		counter++;
		bit_to_check = bit_to_check >> 1;
		//_delay_us(100);	//tx every 0.1ms, 10kHz frequency
		//_delay_us(1100);   // to transmit the signal at 1 kHz frequency
	    //_delay_us(300);
		//_delay_us(300);
	    // _delay_us(400);
		//_delay_us(1);
		//_delay_us(20);
		
		

	}
}

void transmit_32b (unsigned long storage)	// subroutine to transmit SRAM data to PORTA
{
	int counter = 0;
	unsigned long bit_to_check = 0x80000000;  //check the right most bit
	//TIMSK = (1 << TOIE1);		// Timer1 Interrupt Mask Enable
	

	while (counter < 32) 
	{
		if (!(CHECKBIT(storage, bit_to_check)))
		{
			PORTA = 0x00;
      	}
		else
		{
			PORTA = 0x03;
		}

		delay_func(3);
		counter++;
		bit_to_check = bit_to_check >> 1;
		//_delay_us(100);	//tx every 0.1ms, 10kHz frequency
		//_delay_us(1100);   // to transmit the signal at 1 kHz frequency
	    //_delay_us(300);
		//_delay_us(300);
	    // _delay_us(400);
		//_delay_us(1);
		//_delay_us(20);
		
		

	}
}

/*void transmitd(unsigned int data)
{
		int counter = 0, cnt1=0;
		unsigned char bit_to_check=0x80;
		while (counter<8)
		{
			if(!(CHECKBIT(data,bit_to_check)))
			PORTA=0x03;
			counter++;
			bit_to_check = bit_to_check >> 1;
		}

		cnt1 = 0, cnt=0;
		while(cnt1<5)
		{
			cnt1 = cnt;
		}
}*/

void usart_puts( char *str ) 
{ 
    	transmit_ASK(0xFF);
		transmit_ASK(0x80);
     
	 //DDRA=0xFF;
	//PORTA = 0x01;

    while (*str!=0x0a)
	{ 
        transmit_ASK(*str); 
	
        str++; 
		//PORTA ^= 0x01;
    } 
		transmit_ASK(0x00);
	//PORTA ^= 0x01;
}




//-----------------------------------main----------------------------------------



unsigned int ctr;
unsigned char Send[4];

int main(void)
{   

	Timer1_init();
				//sets portD pins as output 
//	DDRB=0xFF;		//sets portB pins as output
	TIMSK |=0x04;	//enable timer1 interrrupt
	DDRA=0xFF;
	PORTA = 0x01;
//	int cnt1=0;
	SREG |= 0x80;	// Enable global interrupt => very important commend!!!
//	unsigned char Rd[32];
Initialize();

    while(1)
    {
			//Get Air Quality Index from PMS sensor
		/*	transmit0(0x42);
			transmit0(0x4d);
			transmit0(0xe2);
			transmit0(0x00);
			transmit0(0x00);
			transmit0(0xFF);
			transmit0(0xFF);
		//	PPM_func(0xAA);
		//	transmit1(0x01);
			for(ctr=0; ctr <32; ctr++)
			{
		//	transmitd(0x7F);	
		//	transmit(0x7F);	
		//	PPM_func(0xFF);
		//	PPM_func(0x7F);
		//	PPM_func(0xFF);
		//	PORTA = 0x0;
		//	
			Rd[ctr] = Rxdata();
		//	transmitd(0x7F);
		//	PPM_func(0x7F);
		
			}*/
		
		//	for(ctr=0; ctr <256; ctr=ctr+4)
						
		//	{
		//	transmit0(0xFF);
		/*	transmit_ASK(0xFF);	//sync byte
			
			transmit_ASK(0x80);	//header byte
		
		    transmit_ASK(0xAA);
			transmit_ASK(0xAB);
			transmit_ASK(0xAC);
			transmit_ASK(0xBB);
			transmit_ASK(0xBC);
			
			transmit_ASK(0x00);*/

       System();
			
		//	}
			
		

  	}	

}

