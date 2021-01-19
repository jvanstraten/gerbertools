/**
 * MIT License
 *
 * Copyright (c) 2021 Jeroen van Straten
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/** \file
 * Handles all the complexity related to aperture macros.
 */

#pragma once

#include <memory>
#include <list>
#include <map>
#include <vector>
#include <sstream>
#include "coord.hpp"
#include "plot.hpp"
#include "aperture.hpp"

namespace aperture_macro {

class Expression;
using ExpressionRef = std::shared_ptr<Expression>;
using ExpressionRefs = std::list<ExpressionRef>;
using Variables = std::map<size_t, double>;

/**
 * Aperture macro expression tree. Can also take the shape of a single token
 * while parsing.
 */
class Expression {
private:

    /**
     * Constructs a string representation of the given expression/token list and
     * context.
     */
    static std::string debug(
        const ExpressionRefs &expr,
        ExpressionRefs::iterator expr_begin,
        ExpressionRefs::iterator expr_end
    );

    /**
     * Reduces the given expression/token list to a single expression tree by
     * applying reduction rules. Yes, it might've been better to just include a
     * proper parser lib and use that. This was easier (by some metrics).
     */
    static ExpressionRef reduce(
        ExpressionRefs &expr,
        ExpressionRefs::iterator expr_begin,
        ExpressionRefs::iterator expr_end
    );

public:
    virtual ~Expression() = default;

    /**
     * Parse the given expression.
     */
    static ExpressionRef parse(std::string expr);

    /**
     * Evaluate this expression with the given set of variables.
     */
    virtual double eval(const Variables &vars) const = 0;

    /**
     * If this is a character token, return the character it represents.
     * Otherwise retun '\0'.
     */
    virtual char get_token() const;

    /**
     * Returns a debug representation of this expression node.
     */
    virtual std::string debug() const = 0;
};

/**
 * Represents a literal.
 */
class LiteralExpression : public Expression {
private:

    /**
     * Value of the literal.
     */
    double value;

public:

    /**
     * Constructs a literal node with the given value.
     */
    explicit LiteralExpression(double value);

    /**
     * Evaluate this expression with the given set of variables.
     */
    double eval(const Variables &vars) const override;

    /**
     * Returns a debug representation of this expression node.
     */
    std::string debug() const override;

};

/**
 * Represents a variable reference.
 */
class VariableExpression : public Expression {
private:

    /**
     * Index of the variable.
     */
    size_t index;

public:

    /**
     * Constructs a variable reference node with the given variable index.
     */
    explicit VariableExpression(size_t index);

    /**
     * Constructs a variable node for the given variable index.
     */
    size_t get_index() const;

    /**
     * Evaluate this expression with the given set of variables.
     */
    double eval(const Variables &vars) const override;

    /**
     * Returns a debug representation of this expression node.
     */
    std::string debug() const override;

};

/**
 * Represents a unary operation, either -x or +x.
 */
class UnaryExpression : public Expression {
private:

    /**
     * The operation, either '-' or '+'.
     */
    char oper;

    /**
     * The operand expression.
     */
    ExpressionRef expr;

public:

    /**
     * Constructs a unary expression node from the given operand character
     * (either '+' or '-') and operand expression.
     */
    UnaryExpression(char oper, const ExpressionRef &expr);

    /**
     * Evaluate this expression with the given set of variables.
     */
    double eval(const Variables &vars) const override;

    /**
     * Returns a debug representation of this expression node.
     */
    std::string debug() const override;

};

/**
 * Represents a unary operation, either x+y, x-y, x*y, or x/y.
 */
class BinaryExpression : public Expression {
private:

    /**
     * The operation, either '-' or '+'.
     */
    char oper;

    /**
     * Left-hand-side operand expression.
     */
    ExpressionRef lhs;

    /**
     * Right-hand-side operand expression.
     */
    ExpressionRef rhs;

public:

    /**
     * Constructs a binary expression node from the given operand character
     * (either '+', '-', 'x', or '/') and operand expressions.
     */
    BinaryExpression(char oper, const ExpressionRef &lhs, const ExpressionRef &rhs);

    /**
     * Evaluate this expression with the given set of variables.
     */
    double eval(const Variables &vars) const override;

    /**
     * Returns a debug representation of this expression node.
     */
    std::string debug() const override;

};

/**
 * Represents an unparsed token.
 */
class Token : public Expression {
private:

    /**
     * The only tokens we need to represent are single characters, so a char is
     * sufficient for identification. These tokens are '(', ')', '+', '-', 'x',
     * and '/'.
     */
    char token;

public:

    /**
     * Constructs a token character. The character must be '(', ')', '+', '-',
     * 'x', or '/'.
     */
    Token(char token);

    /**
     * Evaluate this expression with the given set of variables.
     */
    double eval(const Variables &vars) const override;

    /**
     * Returns the token's character representation.
     */
    char get_token() const override;

    /**
     * Returns a debug representation of this expression node.
     */
    std::string debug() const override;

};

/**
 * Represents an aperture macro command... poorly.
 */
using ApertureMacroCommand = std::vector<ExpressionRef>;

/**
 * Represents an aperture macro (before being instantiated into an aperture).
 */
class ApertureMacro {
private:

    /**
     * The commands that make up the macro.
     */
    std::list<ApertureMacroCommand> cmds;

public:

    /**
     * Parses and appends an aperture macro command.
     */
    void append(const std::string &cmd);

    /**
     * Executes the macro to construct an aperture using the given parameters,
     * reported as a vector of strings. The first string is ignored; it is
     * assumed to be the name of the aperture macro.
     */
    aperture::Ref build(const std::vector<std::string> &csep, const coord::Format &fmt);

};

} // namespace aperture_macro
