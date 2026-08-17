#include "param_parser.h"
#include <cmath>
#include <cctype>
#include <algorithm>
#include <unordered_map>

// ============ Unit Conversion Table ============

namespace {

struct UnitInfo { double toBase; std::string category; };
const std::unordered_map<std::string, UnitInfo> kUnits = {
    // length → m
    {"m",  {1.0,     "length"}},
    {"cm", {0.01,    "length"}},
    {"mm", {0.001,   "length"}},
    {"km", {1000.0,  "length"}},
    // angle → rad
    {"rad", {1.0,        "angle"}},
    {"deg", {M_PI/180.0, "angle"}},
    // voltage → V
    {"V",  {1.0,    "voltage"}},
    {"kV", {1000.0, "voltage"}},
    {"MV", {1e6,    "voltage"}},
    // dimensionless
    {"", {1.0, ""}},
};

double getFactor(const std::string& unit) {
    auto it = kUnits.find(unit);
    return it != kUnits.end() ? it->second.toBase : 1.0;
}

std::string getCategory(const std::string& unit) {
    auto it = kUnits.find(unit);
    return it != kUnits.end() ? it->second.category : "";
}

} // anonymous namespace

// ============ Public API ============

double ParamParser::evaluate(const std::string& expr,
                             const std::map<std::string, double>& vars,
                             const std::string& targetUnit,
                             bool* ok) {
    ParamParser parser(&vars, false, targetUnit);
    parser.tokenize(expr);
    double result = parser.eval();
    if (ok) *ok = parser.m_ok;
    return parser.m_ok ? result : 0.0;
}

std::set<std::string> ParamParser::dependencies(const std::string& expr) {
    std::map<std::string, double> dummy;
    ParamParser parser(&dummy, true, "");
    parser.tokenize(expr);
    parser.eval();
    return parser.m_deps;
}

double ParamParser::convertUnit(double value, const std::string& fromUnit,
                                const std::string& toUnit) {
    if (fromUnit == toUnit) return value;
    std::string catFrom = getCategory(fromUnit);
    std::string catTo   = getCategory(toUnit);
    if (catFrom != catTo || catFrom.empty()) return value; // incompatible or dimensionless
    double factorFrom = getFactor(fromUnit);
    double factorTo   = getFactor(toUnit);
    if (factorFrom == 0.0 || factorTo == 0.0) return value;
    return value * factorFrom / factorTo;
}

// ============ Tokenizer ============

void ParamParser::tokenize(const std::string& expr) {
    m_tokens.clear(); m_pos = 0;
    size_t i = 0;
    while (i < expr.size()) {
        char c = expr[i];
        if (std::isspace(c)) { i++; continue; }
        if (c == '+' || c == '-' || c == '*' || c == '/' || c == '^') {
            m_tokens.push_back(Token(Token::OP, std::string(1, c)));
            i++;
        } else if (c == '(') { m_tokens.push_back(Token(Token::LPAREN)); i++; }
        else if (c == ')') { m_tokens.push_back(Token(Token::RPAREN)); i++; }
        else if (c == ',') { m_tokens.push_back(Token(Token::COMMA)); i++; }
        else if (c == '[') {
            // unit specifier: [...]
            i++;
            size_t start = i;
            while (i < expr.size() && expr[i] != ']') i++;
            std::string unit = expr.substr(start, i - start);
            if (i < expr.size()) i++; // skip ']'
            if (!m_tokens.empty()) m_tokens.back().unit = unit;
        } else if (std::isdigit(c) || c == '.') {
            size_t end; double v = std::stod(expr.substr(i), &end);
            m_tokens.push_back(Token(v));
            i += end;
        } else if (std::isalpha(c) || c == '_') {
            size_t start = i;
            while (i < expr.size() && (std::isalnum(expr[i]) || expr[i] == '_')) i++;
            std::string name = expr.substr(start, i - start);
            m_tokens.push_back(Token(Token::ID, name));
        } else { i++; }
    }
    m_tokens.push_back(Token(Token::END));
}

// ============ Recursive Descent Parser ============

double ParamParser::applyUnit(double value) {
    if (cur().unit.empty() || m_targetUnit.empty()) return value;
    return convertUnit(value, cur().unit, m_targetUnit);
}

double ParamParser::eval() {
    m_pos = 0; m_ok = true;
    double r = parseExpr();
    if (cur().type != Token::END) m_ok = false;
    return r;
}

double ParamParser::parseExpr() {
    double left = parseTerm();
    while (m_ok && cur().type == Token::OP &&
           (cur().name == "+" || cur().name == "-")) {
        std::string op = cur().name; advance();
        double right = parseTerm();
        left = (op == "+") ? left + right : left - right;
    }
    return left;
}

double ParamParser::parseTerm() {
    double left = parsePower();
    while (m_ok && cur().type == Token::OP &&
           (cur().name == "*" || cur().name == "/")) {
        std::string op = cur().name; advance();
        double right = parsePower();
        if (op == "/" && std::abs(right) < 1e-15) { m_ok = false; return 0; }
        left = (op == "*") ? left * right : left / right;
    }
    return left;
}

double ParamParser::parsePower() {
    double left = parseUnary();
    while (m_ok && cur().type == Token::OP && cur().name == "^") {
        advance();
        double right = parseUnary();
        left = std::pow(left, right);
    }
    return left;
}

double ParamParser::parseUnary() {
    if (cur().type == Token::OP && cur().name == "-") {
        advance();
        return -parseAtom();
    }
    return parseAtom();
}

double ParamParser::parseAtom() {
    if (cur().type == Token::NUM) {
        double v = cur().val;
        v = applyUnit(v);  // convert 5[cm] → 0.05 (if targetUnit is "m")
        advance(); return v;
    }
    if (cur().type == Token::LPAREN) {
        advance();
        double v = parseExpr();
        if (cur().type != Token::RPAREN) { m_ok = false; return 0; }
        advance();
        v = applyUnit(v); // (expr)[cm]
        return v;
    }
    if (cur().type == Token::ID) {
        std::string name = cur().name;
        if (name == "pi") { advance(); return M_PI; }
        if (name == "e")  { advance(); return M_E; }
        advance();
        if (cur().type == Token::LPAREN) {
            advance();
            double arg = parseExpr();
            if (name == "sqrt")  arg = std::sqrt(arg);
            else if (name == "abs")  arg = std::abs(arg);
            else if (name == "sin")  arg = std::sin(arg);
            else if (name == "cos")  arg = std::cos(arg);
            else if (name == "tan")  arg = std::tan(arg);
            else if (name == "exp")  arg = std::exp(arg);
            else if (name == "log")  arg = std::log(arg);
            else if (name == "min" || name == "max") {
                if (cur().type != Token::COMMA) { m_ok = false; return 0; }
                advance();
                double arg2 = parseExpr();
                arg = (name == "min") ? std::min(arg, arg2) : std::max(arg, arg2);
            } else if (name == "pow") {
                if (cur().type != Token::COMMA) { m_ok = false; return 0; }
                advance();
                double arg2 = parseExpr();
                arg = std::pow(arg, arg2);
            } else { m_ok = false; return 0; }
            if (cur().type != Token::RPAREN) { m_ok = false; return 0; }
            advance();
            return arg;
        }
        // 变量引用
        if (m_collectDeps) { m_deps.insert(name); return 0.0; }
        if (m_vars) {
            auto it = m_vars->find(name);
            if (it != m_vars->end()) return it->second;
        }
        m_ok = false; return 0;
    }
    m_ok = false; return 0;
}
