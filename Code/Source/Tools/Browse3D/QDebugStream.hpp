#ifndef Q_DEBUGSTREAM_H
#define Q_DEBUGSTREAM_H

#include <QtGlobal>
#include <iostream>
#include <streambuf>
#include <string>

/*
    A convenience class for wrapping the standard iostreams into
    objects, which will send the output to qdebug.

    It is faster, of course, to write directly to qDebug, but if you
    have existing code that uses iostream and want a cheap way of redirecting,
    this is it.

    Usage:
    QApplication app(argc, argv);

    // these redirect both cout/cerr
    QDebugStream qdsOut(std::cout);
    QDebugStream qdsErr(std::cerr);

    // now start using cout and cerr normally
    std::cerr << "Oops";   // this goes to your QDebug handler

    http://www.qtforum.org/article/678/redirecting-cout-cerr-to-qdebug.html
*/
class QDebugStream : public std::basic_streambuf<char>
{
public:
  QDebugStream(std::ostream &stream) : m_stream(stream)
  {
    m_old_buf = stream.rdbuf();
    stream.rdbuf(this);
  }

  ~QDebugStream()
  {
    // output anything that is left
    if (!m_string.empty())
      qDebug("%s", m_string.c_str());

    m_stream.rdbuf(m_old_buf);
  }

protected:
  virtual int_type overflow(int_type v)
  {
    if (v == '\n')
    {
      qDebug("%s", m_string.c_str());
      m_string.clear();
    }
    else
      m_string.push_back(v);

    return v;
  }

  virtual std::streamsize xsputn(const char *p, std::streamsize n)
  {
    m_string.append(p, p + n);

    std::size_t pos = 0;
    while (pos != std::string::npos)
    {
      pos = m_string.find('\n');
      if (pos != std::string::npos)
      {
        std::string tmp(m_string.begin(), m_string.begin() + pos);
        qDebug("%s", tmp.c_str());
        m_string.erase(m_string.begin(), m_string.begin() + pos + 1);
      }
    }

    return n;
  }

private:
  std::ostream &m_stream;
  std::streambuf *m_old_buf;
  std::string m_string;
};

#endif
