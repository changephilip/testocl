#ifndef _C_STRING_CUH_
#define _C_STRING_CUH_


/* Provide an uniform c string library.
   All following functions can be found by internet.
*/

bool c_strcmp(const char *a, const char *b, unsigned l) {
    for (int i = 0; i< l; i++)
        if (a[i] != b[i]) return false;
    return true;
}

bool c_strncmp(const char* s1, const char* s2, size_t n) {
    while(n--)
        if(*s1++!=*s2++)
            return *(unsigned char*)(s1 - 1) - *(unsigned char*)(s2 - 1);
    return 0;
}

int c_strspn(const char *s, const char *accept) {
    const char *p, *a;
    size_t count = 0;

    for (p = s; *p != '\0'; ++p) {
        for (a = accept; *a != '\0'; ++a)
            if (*p == *a)
                break;
        if (*a == '\0')
            return count;
        else
            ++count;
    }

    return count;
}

int c_strlen(const char * str) {
    const char *s;
    for (s = str; *s; ++s) {}
    return(s - str);
}

char *c_strchr(char *s, int c) {
    while(*s != '\0' && *s != c) { ++s; }
    return *s == c ? s : NULL;
}

char *c_strpbrk(const char *cs, const char *ct) {
    const char *sc1, *sc2;
    for( sc1 = cs; *sc1 != '\0'; ++sc1) {
        for( sc2 = ct; *sc2 != '\0'; ++sc2) {
            if (*sc1 == *sc2) {
                return (char *) sc1;
            }
        }
    }
    return NULL;
}

char *c_strtok_r(char *s, const char *delim, char **save_ptr) {
    char *token;

    if (s == NULL) s = *save_ptr;

    s += c_strspn(s, delim);
    if (*s == '\0')
        return NULL;

    token = s;
    s = c_strpbrk(token, delim);
    if (s == NULL)
        *save_ptr = c_strchr(token, '\0');
    else {
        *s = '\0';
        *save_ptr = s + 1;
    }

    return token;
}

void c_strcpy(char *dest, const char *source) {
    int i = 0;
    while ((dest[i] = source[i]) != '\0') { i++; }
}

char *c_strncpy(char *dest, const char *src, size_t n) {
    char *ret = dest;
    do {
        if (!n--)
            return ret;
    } while (*dest++ = *src++);

    while (n--)
        *dest++ = 0;

    return ret;
}

void *c_memset(void *s, int c, size_t n) {
    unsigned char* p = static_cast<unsigned char *>(s);
    while(n--)
        *p++ = (unsigned char)c;
    return s;
}

bool is_digitstr(char *str) {
    return (c_strspn(str, "0123456789") == c_strlen(str));
}

int c_atoi(char *s) {
    int acum = 0;
    while((*s >= '0')&&(*s <= '9')) {
        acum = acum * 10;
        acum = acum + (*s - 48);
        s++;
    }
    return acum;
}

#endif //_CUDA_STRING_CUH_
