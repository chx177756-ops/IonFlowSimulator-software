#pragma once
#include <string>
#include <map>
#include <set>
#include <vector>

// 递归下降表达式解析器, 不依赖Qt/OCC
class ParamParser {
public:
    // 求值表达式, vars提供参数名→值映射. targetUnit为期望输出单位.
    static double evaluate(const std::string& expr,
                           const std::map<std::string, double>& vars,
                           const std::string& targetUnit = "",
                           bool* ok = nullptr);
    // 提取表达式中引用的参数名列表
    static std::set<std::string> dependencies(const std::string& expr);
    // 单位转换
    static double convertUnit(double value, const std::string& fromUnit,
                              const std::string& toUnit);

private:
    struct Token {
        enum Type { NUM, ID, OP, LPAREN, RPAREN, COMMA, END };
        Type type; double val; std::string name; std::string unit;
        Token(Type t) : type(t), val(0) {}
        Token(double v) : type(NUM), val(v) {}
        Token(Type t, const std::string& s) : type(t), name(s) {}
    };

    const std::map<std::string, double>* m_vars = nullptr;
    std::vector<Token> m_tokens;
    size_t m_pos = 0;
    bool m_ok = true;
    bool m_collectDeps = false;
    std::set<std::string> m_deps;
    std::string m_targetUnit;

    ParamParser(const std::map<std::string, double>* vars, bool collectDeps,
                const std::string& targetUnit)
        : m_vars(vars), m_collectDeps(collectDeps), m_targetUnit(targetUnit) {}

    void tokenize(const std::string& expr);
    Token& cur() { return m_tokens[m_pos]; }
    void advance() { if (m_pos < m_tokens.size()) m_pos++; }

    double parseExpr();
    double parseTerm();
    double parsePower();
    double parseUnary();
    double parseAtom();   // returns value in targetUnit

    double eval();
    double applyUnit(double value);  // convert value from cur().unit to targetUnit
};
