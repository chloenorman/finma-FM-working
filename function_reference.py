"""
FINANCE & GENERAL USE - FUNCTION AND METHOD REFERENCE
======================================================
This module provides a comprehensive summary of important functions and methods
available in this repository for financial calculations and general Python use.

Repository Structure:
- ytm/: Yield-to-Maturity and yield curve analysis
- intro_coding/: Pandas and portfolio analysis examples
"""

# ============================================================================
# YIELD CURVE & INTEREST RATE CLASSES & METHODS
# ============================================================================
# Located in: ytm/curve_classes_and_functions.py

class ZeroCurve:
    """
    Stores zero rates and discount factors for a given set of maturities.
    
    Methods:
    --------
    add_zero_rate(maturity, zero_rate)
        - Adds a zero rate to the curve
    
    add_discount_factor(maturity, discount_factor)
        - Adds a discount factor to the curve
    
    get_AtMat(maturity)
        - Returns the accumulated value (e^(r*T)) at a given maturity
        - Uses exponential interpolation if maturity not in curve
    
    get_discount_factor(maturity)
        - Returns the discount factor at a given maturity
        - Uses exponential interpolation if not explicitly stored
    
    get_zero_rate(maturity)
        - Returns the zero rate (spot rate) for a given maturity
    
    get_zero_curve()
        - Returns tuple of (maturities, discount_factors)
    
    npv(cash_flows)
        - Calculates Net Present Value of cash flows using the yield curve
        - Parameters: CashFlows object
        - Returns: float (NPV)
    """
    pass


class YieldCurve(ZeroCurve):
    """
    Extends ZeroCurve to build a yield curve through bootstrapping.
    Used for extracting spot rates from bond prices.
    
    Methods:
    --------
    set_constituent_portfolio(portfolio)
        - Sets the portfolio of bills and bonds used for bootstrapping
    
    bootstrap()
        - Bootstraps the yield curve from a portfolio of bills and bonds
        - Extracts zero rates from market prices
        - Bills are priced first, then bonds are used for longer maturities
    """
    pass


def exp_interp(xs, ys, x):
    """
    Exponential interpolation for continuously compounded rates.
    
    Parameters:
    -----------
    xs : list or np.array
        Vector of x values (sorted)
    ys : list or np.array
        Vector of y values (accumulated factors)
    x : float
        The x value to interpolate
    
    Returns:
    --------
    float
        Interpolated y value using continuous compounding
    
    Use Case:
    ---------
    Used internally by ZeroCurve to interpolate discount factors and rates
    for maturities not explicitly in the curve.
    """
    pass


# ============================================================================
# FIXED INCOME INSTRUMENT CLASSES
# ============================================================================
# Located in: ytm/instrument_classes.py

class CashFlows:
    """
    Base class for managing cash flows with associated maturities.
    
    Methods:
    --------
    add_cash_flow(maturity, amount)
        - Adds a cash flow at a specific maturity
    
    get_cash_flow(maturity)
        - Returns the cash flow at a given maturity
        - Returns None if maturity not found
    
    get_maturities()
        - Returns list of all maturity dates
    
    get_amounts()
        - Returns list of all cash flow amounts
    
    get_cash_flows()
        - Returns list of (maturity, amount) tuples
    """
    pass


class Bank_bill(CashFlows):
    """
    Short-term debt instrument (typically < 1 year).
    
    Parameters (Constructor):
    -------------------------
    face_value : float, default=100
        Par value of the bill
    maturity : float, default=0.25
        Time to maturity in years
    ytm : float, default=0.00
        Yield to maturity
    price : float, default=100
        Current market price
    
    Methods:
    --------
    set_face_value(face_value) / get_face_value()
    set_maturity(maturity) / get_maturity()
    set_ytm(ytm) / get_ytm()
        - Sets YTM and updates price automatically
    set_price(price)
        - Sets price and updates YTM automatically
    get_price()
    set_cash_flows()
        - Generates cash flows: [-price at t=0, +face_value at maturity]
    
    Pricing Formula:
    ----------------
    Price = FaceValue / (1 + YTM * maturity)
    """
    pass


class Bond(CashFlows):
    """
    Fixed-rate coupon bond.
    
    Parameters (Constructor):
    -------------------------
    face_value : float, default=100
        Par value of the bond
    maturity : float, default=3
        Time to maturity in years
    coupon : float, default=0
        Annual coupon rate (as decimal, e.g., 0.05 for 5%)
    frequency : int, default=4
        Coupon payment frequency per year (e.g., 4 for quarterly)
    ytm : float, default=0
        Yield to maturity
    price : float, default=100
        Current market price
    
    Methods:
    --------
    set_face_value(face_value) / get_face_value()
    set_maturity(maturity) / get_maturity()
    set_coupon(coupon) / get_coupon_rate()
    set_frequency(frequency)
    set_ytm(ytm) / get_ytm()
        - Sets YTM and updates price automatically
    get_price()
    set_cash_flows()
        - Generates cash flows for all coupon and maturity payments
    
    Pricing Formula:
    ----------------
    Bond Price = Sum of PV(coupon payments) + PV(face value)
    = [C/freq * (1 - (1+y/freq)^(-n*freq)) / (y/freq)] + FV / (1+y/freq)^(n*freq)
    where:
        C = annual coupon rate
        y = YTM
        freq = payment frequency
        n = years to maturity
    """
    pass


class Portfolio(CashFlows):
    """
    Collection of bonds and bank bills.
    
    Methods:
    --------
    add_bond(bond)
        - Adds a Bond to the portfolio
    
    add_bank_bill(bank_bill)
        - Adds a Bank_bill to the portfolio
    
    get_bonds()
        - Returns list of all bonds
    
    get_bank_bills()
        - Returns list of all bank bills
    
    set_cash_flows()
        - Aggregates all cash flows from all instruments
    
    Use Case:
    ---------
    Used in YieldCurve.bootstrap() to extract zero rates from market prices
    """
    pass


# ============================================================================
# COMMON FINANCIAL CALCULATIONS & PATTERNS
# ============================================================================

def calculate_npv(rate, cash_flows):
    """
    Calculate Net Present Value of a series of cash flows.
    
    Usage:
    ------
    >>> from ytm.curve_classes_and_functions import ZeroCurve
    >>> curve = ZeroCurve()
    >>> curve.add_zero_rate(1, 0.05)
    >>> curve.add_zero_rate(2, 0.06)
    >>> cf = CashFlows()
    >>> cf.add_cash_flow(1, 100)
    >>> cf.add_cash_flow(2, 100)
    >>> npv = curve.npv(cf)
    """
    pass


def bootstrap_yield_curve(bills, bonds):
    """
    Build a zero-coupon yield curve from market prices of bills and bonds.
    
    Parameters:
    -----------
    bills : list of Bank_bill
        Sorted by maturity (ascending)
    bonds : list of Bond
        Sorted by maturity (ascending)
        Each bond should introduce one new maturity beyond previous instruments
    
    Returns:
    --------
    YieldCurve object with zero rates for all maturities
    
    Process:
    --------
    1. Use bill prices to extract zero rates for short maturities
    2. Use bond prices to extract zero rates for longer maturities
    3. Use previously calculated discount factors to back out new rates
    """
    pass


# ============================================================================
# PANDAS & DATA ANALYSIS (from intro_coding notebooks)
# ============================================================================

# Common pandas operations used in this repository:

def import_portfolio_data():
    """
    Import and analyze portfolio holdings and transactions.
    
    Files Used:
    -----------
    - holdings.csv: Current portfolio composition
    - transactions.csv: Historical portfolio transactions
    - company_info.csv: Corporate information
    - stock_prices_timeseries.csv: Historical price data
    
    Typical Operations:
    -------------------
    1. Read CSV files into pandas DataFrames
    2. Calculate daily returns from price data
    3. Analyze portfolio weights and composition
    4. Compute portfolio summary statistics (return, volatility)
    5. Generate reports (portfolio_summary.csv)
    """
    pass


# ============================================================================
# QUICK REFERENCE: KEY FORMULAS
# ============================================================================

"""
DISCOUNT FACTOR (Zero-Coupon Bond):
    DF(T) = e^(-r*T)
    where r = spot rate (zero rate), T = time to maturity

PRESENT VALUE:
    PV = CF / DF(T) = CF * e^(r*T)

BOND PRICE (Coupon Bond):
    P = Sum[C/freq / (1 + y/freq)^(t*freq) for t in periods] + FV / (1 + y/freq)^(n*freq)

BANK BILL PRICE:
    P = FV / (1 + y*T)
    where y = simple annual rate, T = fraction of year

CONTINUOUSLY COMPOUNDED RATES:
    r_cont = ln(AccumulationFactor)
    AccumulationFactor = e^(r_cont * T)

INTERPOLATION (Exponential):
    For rates: Use ln(y) as linear interpolation, then exponentiate
    For accumulated factors: y = y0 * exp(rate * (x - x0))
"""


# ============================================================================
# WORKFLOW EXAMPLE: BOND PRICING & YIELD CURVE
# ============================================================================

"""
1. CREATE INSTRUMENTS:
   bill = Bank_bill(face_value=100, maturity=0.25, price=98.75)
   bond = Bond(face_value=100, maturity=2, coupon=0.05, frequency=2, price=102.5)

2. BUILD PORTFOLIO:
   portfolio = Portfolio()
   portfolio.add_bank_bill(bill)
   portfolio.add_bond(bond)

3. BOOTSTRAP YIELD CURVE:
   curve = YieldCurve()
   curve.set_constituent_portfolio(portfolio)
   curve.bootstrap()

4. USE CURVE FOR VALUATION:
   new_bond = Bond(face_value=100, maturity=1, coupon=0.04, frequency=2)
   new_bond.set_cash_flows()
   npv = curve.npv(new_bond.get_cash_flows())

5. INTERPOLATE FOR OTHER MATURITIES:
   df_18m = curve.get_discount_factor(1.5)
   rate_18m = curve.get_zero_rate(1.5)
"""

