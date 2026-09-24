
#   include <string>
#   include <sstream>

#   pragma once

#   ifndef __STRING_UTILS__
#   define __STRING_UTILS__

    /*
    --------------------------------------------------------
     * SEPARATOR: return the OS-specific path separator
    --------------------------------------------------------
     */

    inline char separator (
        )
    {
#   ifdef _WIN32
        return '\\';
#   else
        return '/';
#   endif
    }

    /*
    --------------------------------------------------------
     * PATH-JOIN: joint two path names with a separator
    --------------------------------------------------------
     */

    std::string path_join (
        std::string const& _root,
        std::string const& _stem
        )
    {
        if(!_root.empty() && !_stem.empty() &&
            _root. back() != separator() &&
            _stem.front() != separator())
        {
            return _root + separator() + _stem;
        }
        else
        {
            return _root + _stem;
        }
    }

    /*
    --------------------------------------------------------
     * FILE-PART: split a file name into path/name.fext
    --------------------------------------------------------
     */

    void file_part (
        std::string const& _fsrc,
        std::string      & _path,
        std::string      & _name,
        std::string      & _fext
        )
    {
        typename std::string::size_type
            _spos = _fsrc.find_last_of("\\/");

        typename std::string::size_type
            _dpos = _fsrc.find_last_of("."  );

        typename std::string::const_iterator
            _pos0, _pos1, _pos2,
            _pos3, _pos4, _pos5;

        if (_spos != std::string::npos )
        {
            _pos0 = _fsrc.begin();
            _pos1 = _fsrc.begin()+_spos-0;
            _pos2 = _fsrc.begin()+_spos+1;
        }
        else
        {
            _pos0 = _fsrc.begin();
            _pos1 = _fsrc.begin();
            _pos2 = _fsrc.begin();
        }

        if (_dpos != std::string::npos &&
           (_spos == std::string::npos ||
            _dpos >= _spos) )
        {
            _pos3 = _fsrc.begin()+_dpos-0;
            _pos4 = _fsrc.begin()+_dpos+1;
            _pos5 = _fsrc.end  ();
        }
        else
        {
            _pos3 = _fsrc.end  ();
            _pos4 = _fsrc.end  ();
            _pos5 = _fsrc.end  ();
        }

        _path = std::string(_pos0, _pos1);
        _name = std::string(_pos2, _pos3);
        _fext = std::string(_pos4, _pos5);
    }

    /*
    --------------------------------------------------------
     * CF-CONVENTIONS: the Conventions attribute for output,
     * keeping the entries (including any CF version) from
     * the input and adding "CF-1.8" and "MPAS" if missing
    --------------------------------------------------------
     */

    inline std::string cf_conventions (
        std::string const& _in  // input Conventions, may be empty
        )
    {
        std::string _list = _in, _item, _cf, _rest;
        bool _mpas = false;

        // entries are separated by blanks or commas; also drop the
        // trailing null from reading the attribute
        for (auto &_char : _list)
        {
            if (_char == ',' || _char == '\0') _char = ' ';
        }

        std::istringstream _stream(_list);
        while (_stream >> _item)
        {
            if (_item.compare(0, 3, "CF-") == 0)
            {
                if (_cf.empty()) _cf = _item;
                continue;
            }
            if (_item == "MPAS") _mpas = true;
            _rest += " " + _item;
        }

        if (_cf.empty()) _cf = "CF-1.8";
        if (!_mpas) _rest += " MPAS";

        return _cf + _rest;
    }

#   endif   //__STRING_UTILS__
