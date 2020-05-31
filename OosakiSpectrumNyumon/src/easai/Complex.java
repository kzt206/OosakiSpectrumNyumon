//Complex.java  -- An implementation of a complex number class.

//Copyright (c) 2013 easai

//Author: easai 
//Website: easai.web.fc2.com/homepage/javadev/Trig/
//Created: Wed Dec 18 06:17:27 2013
//Keywords: math, mathematics, complex numbers, imaginary numbers, complex analysis

//Commentary:
//This is an implementation of a complex number class.
//
//Disclaimer:
//There is no way to overload operators in Java, the expressions are rather clumsy, z.add(z0) etc.
//

//Code:

//package com.github.easai.math;
package easai;

/**
 * The <tt>Complex</tt> class represents a complex number class.
 * 
 * @author easai
 */
public class Complex {
	/**
	 * The real part of the complex number.
	 */
	private double real = 0;
	/**
	 * The imaginary part of the complex number.
	 */
	private double imaginary = 0;

	/**
	 * The imaginary number i.
	 */
	public static final Complex I = new Complex(0, 1);

	/**
	 * Constructs a complex number.
	 * 
	 * @param real the real part of the complex number
	 * @param imaginary the imaginary part of the complex number
	 */
	public Complex(double real, double imaginary) {
		this.real = real;
		this.imaginary = imaginary;
	}

	/**
	 * Returns the real part of the complex number.
	 * 
	 * @return the real part of the complex number.
	 */
	public double Re() {
		return real;
	}

	/**
	 * Returns the imaginary part of the complex number.
	 * 
	 * @return the imaginary part of the complex number.
	 */
	public double Im() {
		return imaginary;
	}

	/**
	 * Adds the complex number.
	 * 
	 * @param c
	 *            the complex number that will be added
	 * @return the sum (a new instance of Complex class)
	 */
	public Complex add(Complex c) {
		return new Complex(real + c.real, imaginary + c.imaginary);
	}

	/**
	 * Subtracts the complex number.
	 * 
	 * @param c
	 *            the complex number that will be subtracted
	 * @return the difference (a new instance of Complex class)
	 */
	public Complex subtract(Complex c) {
		return new Complex(real - c.real, imaginary - c.imaginary);
	}

	/**
	 * Multiplies the complex number.
	 * 
	 * @param c the complex number that will be multiplied
	 * @return the product (a new instance of Complex class)
	 */
	public Complex multiply(Complex c) {
		// (x+iy)*(xx+iyy)=x*xx+x*iyy+iy*xx-y*yy
		// for debugging purposes only
		// System.out.println(x+"*"+c.x+"-"+y+"*"+c.y+" = "+(x*c.x-y*c.y));
		return new Complex(real * c.real - imaginary * c.imaginary, real * c.imaginary + imaginary * c.real);
	}

	/**
	 * Returns (-1)*z.
	 * 
	 * @return (-1)*z (a new instance of Complex class)
	 */
	public Complex minus() {
		return multiply(new Complex(-1, 0));
	}

	/**
	 * Returns true if the complex number has the same values.
	 * 	
	 * @param c the complex number that will be compared with  
	 * @return returns true if the complex number has the same values.
	 */
	public boolean equals(Complex c) {
		return (real == c.real) && (imaginary == c.imaginary);
	}

	/**
	 * Returns a string representation of the complex number.
	 * 
	 * @return a string representation of the complex number.
	 */
	public String toString() {
		return real + " + " + imaginary + "*I";
	}

	/**
	 * Returns the absolute values of the complex number, i.e. |z|.
	 * 
	 * @return the absolute values of the complex number, i.e. |z|.
	 */
	public double abs() {
		return Math.sqrt(real * real + imaginary * imaginary);
	}

	/**
	 * Returns the value of e to the complex number, i.e. e^z.
	 * 
	 * @param z the exponent
	 * @return the value of e to the complex number, i.e. e^z.
	 */
	public static Complex exp(Complex z) {
		// z=x+iy
		// e^x*e^iy
		// e^iy=cosy+i*siny
		// Re e^z=e^x*cosy
		// Im e^z=e^x*siny
		double ex = Math.pow(Math.E, z.real);
		return new Complex(ex * Math.cos(z.imaginary), ex * Math.sin(z.imaginary));
	}

	/**
	 * Returns the sine value of the complex number, i.e. sin(z).
	 * 
	 * @param z the 'angle' value (a complex number)
	 * @return the sine value of the complex number, i.e. sin(z).
	 */
	public static Complex sin(Complex z) {
		Complex zi = z.multiply(Complex.I);
		Complex diff = exp(zi.minus()).subtract(exp(zi));
		return diff.multiply(new Complex(.5, 0)).multiply(Complex.I);
	}

	/**
	 * Returns cosine of the complex number, i.e. cos(z).
	 * 
	 * @param z the 'angle' value (a complex number)
	 * @return the cosine value of the complex number, i.e. sin(z).
	 */
	public static Complex cos(Complex z) {
		Complex zi = z.multiply(Complex.I);
		Complex sum = exp(zi).add(exp(zi.minus()));
		return sum.multiply(new Complex(.5, 0));
	}

	/**
	 * If successfully complied, it should output the following.
	 * 
	 * <pre>
	 * java Complex
	 * z = 0.0 + 1.0*I
	 * sin(z) = 0.0 + 1.1752011936438014*I
	 * cos(z) = 1.5430806348152437 + 0.0*I
	 * e^z = 0.5403023058681398 + 0.8414709848078965*I
	 * </pre>
	 * 
	 * @param args not used
	 */
	public static void main(String args[]) {
		Complex z = new Complex(0, 1);
		Complex sin = Complex.sin(z);
		Complex cos = Complex.cos(z);
		Complex exp = Complex.exp(z);
		System.out.println("z = " + z.toString());
		System.out.println("sin(z) = " + sin.toString());
		System.out.println("cos(z) = " + cos.toString());
		System.out.println("e^z = " + exp.toString());
	}
}

// Complex.java ends here